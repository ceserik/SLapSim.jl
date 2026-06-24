
const _N_ELLIPSE_PTS = 41
const _ELLIPSE_ANGLES = range(0, 2π, length=_N_ELLIPSE_PTS)

mutable struct  Tire{F1,F2,F3,F4,F5}
    radius::carParameter{carVar}
    width::carParameter{carVar}
    inertia::carParameter{carVar}
    mass::carParameter{carVar}
    velocity::carParameter{Vector{carVar}}
    angularFrequency::carParameter{carVar}
    forces::carParameter{Vector{carVar}}
    slipAngle::carParameter{carVar}
    slipRatio::carParameter{carVar}
    compute::F1
    tireConstraints::F2
    setVelocity::F3
    maxSLipAngle::carParameter{carVar}
    scalingForce::carParameter{carVar}
    frictionCoefficient::carParameter{carVar}
    rollingResistance::carParameter{carVar}
    brakingForce::carParameter{carVar}
    setupObservables::F4
    updateObservables::F5
end
Base.show(io::IO, ::MIME"text/plain", obj::Tire) = prettyPrintComponent(io, obj)

function plotTire(tire::Tire; Fz::Float64=1000.0, max_slip_deg::Float64=15.0, label::String="", ax=nothing)
    αs = range(-max_slip_deg, max_slip_deg, length=200) .* (π / 180)
    Fy = similar(αs)

    old_forces = copy(tire.forces.value)
    old_vel = copy(tire.velocity.value)

    for (i, α) in enumerate(αs)
        tire.forces.value[3] = Fz
        vx = 10.0
        tire.velocity.value .= [vx, -vx * tan(α), 0.0]
        tire.compute(0.0)
        Fy[i] = tire.forces.value[2]
    end

    tire.forces.value .= old_forces
    tire.velocity.value .= old_vel

    αs_deg = αs .* (180 / π)
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel="Slip angle [deg]", ylabel="Lateral force [N]", title="Tire curve (Fz=$(Fz)N)")
        lines!(ax, αs_deg, Fy; linewidth=2, label=label)
        display(GLMakie.Screen(), fig)
        return fig
    else
        lines!(ax, αs_deg, Fy; linewidth=2, label=label)
        return ax
    end
end

function createR20lin(motor::Motor,gearbox::Gearbox)
    radius = carParameter{carVar}(0.205, "tire radius", "m")
    width = carParameter{carVar}(0.3, "tire width, wrong", "m")
    inertia = carParameter{carVar}(0.3, "tire width, wrong", "m")
    mass = carParameter{carVar}(1.0, "tire mass,wrong", "kg")
    velocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "velocity", "m/s")
    angularFrequency = carParameter{carVar}(0.0, "angular velocity", "rad/s")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Tire Force", "N")
    slipAngle = carParameter{carVar}(0.0, "Slip angle", "rad")
    slipRatio = carParameter{carVar}(0.0, "slip ratio", "-")
    maxSlipAngle = carParameter{carVar}(5/180*pi, "max slip angle", "rad")
    scalingForce = carParameter{carVar}(0.0, "max longitudinal force", "N")
    scalingForce.value = motor.torqueSpeedFunction(0.0) * gearbox.ratio.value / radius.value
    frictionCoefficient = carParameter{carVar}(1.0, "Friction Coefficient", "-", :sensitivity)
    rollingResistance = carParameter{carVar}(0.010, "Rolling resistance coeff", "-", :sensitivity)
    brakingForce = carParameter{carVar}(0.0, "braking force", "N")

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
        #static (non-dynamic) tire: angular frequency follows longitudinal velocity.
        #dynamic tire variants should override this and treat angularFrequency as a state.
        angularFrequency.value = velocityIn[1] / radius.value
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])

        if isnothing(optiModel)
            forces.value[2] = 1*slipAngle.value/maxSlipAngle.value * forces.value[3] * frictionCoefficient.value
            forces.value[1] = inTorque/radius.value + brakingForce.value
        else
            sf = scalingForce.value
            fy_expr = (slipAngle.value/maxSlipAngle.value * forces.value[3] * frictionCoefficient.value) / sf
            fx_expr = (inTorque/radius.value + brakingForce.value) / sf
            fy_raw = @variable(optiModel, lower_bound=-2.0, upper_bound=2.0, start=0.0)
            fx_raw = @variable(optiModel, lower_bound=-2.0, upper_bound=2.0, start=0.0)
            @constraint(optiModel, fy_raw == fy_expr)
            @constraint(optiModel, fx_raw == fx_expr)
            forces.value[2] = fy_raw * sf
            forces.value[1] = fx_raw * sf
        end
    end

    function tireConstraints(model=nothing)
        Fx_s = forces.value[1]/scalingForce.value
        Fy_s = forces.value[2]/scalingForce.value
        Fz_s = frictionCoefficient.value * forces.value[3]/scalingForce.value
        if isnothing(model)
            lessContraint(Fx_s^2 + Fy_s^2, Fz_s^2, nothing)
        else
            @constraint(model, Fx_s^2 + Fy_s^2  <= Fz_s^2)
        end
    end

    function setupObservables(ax::Axis)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        rect_obs = Observable(dummy)
        poly!(ax, rect_obs; color=(:black, 0.3), strokecolor=(:gray30, 0.5), strokewidth=1)
        ellipse_obs = Observable([Point2f(0, 0) for _ in _ELLIPSE_ANGLES])
        lines!(ax, ellipse_obs; color=(:gray, 0.5), linewidth=1, linestyle=:dash)
        force_pos_obs = Observable([Point2f(0, 0)])
        force_dir_obs = Observable([Point2f(0, 0)])
        arrows2d!(ax, force_pos_obs, force_dir_obs; color=:green, shaftwidth=4, tipwidth=8, tiplength=8)
        return (rect=rect_obs, ellipse=ellipse_obs, force_pos=force_pos_obs, force_dir=force_dir_obs)
    end

    function updateObservables(obs::NamedTuple, wx::Float64, wy::Float64, θ::Float64, Fz_static::Float64)
        r = radius.value
        s = 4r / Fz_static
        R = _rotmat2d(θ)
        obs.rect[] = _rect_points(wx, wy, 2r, width.value * 0.8, θ)
        Fmax_s = abs(forces.value[3]) * frictionCoefficient.value * s
        obs.ellipse[] = [Point2f(R * [Fmax_s * cos(a), Fmax_s * sin(a)] + [wx, wy]) for a in _ELLIPSE_ANGLES]
        f_global = R * [forces.value[1], forces.value[2]] .* s
        obs.force_pos[] = [Point2f(wx, wy)]
        obs.force_dir[] = [Point2f(f_global[1], f_global[2])]
    end

    tire = Tire(
        radius,
        width,
        inertia,
        mass,
        velocity,
        angularFrequency,
        forces,
        slipAngle,
        slipRatio,
        compute,
        tireConstraints,
        setVelocity,
        maxSlipAngle,
        scalingForce,
        frictionCoefficient,
        rollingResistance,
        brakingForce,
        setupObservables,
        updateObservables,
    )
end



function createR20_pacejka(motor::Motor,gearbox::Gearbox)
    radius = carParameter{carVar}(0.205, "tire radius", "m")
    width = carParameter{carVar}(0.3, "tire width, wrong", "m")
    inertia = carParameter{carVar}(0.3, "tire inertia, wrong", "m")
    mass = carParameter{carVar}(1.0, "tire mass,wrong", "kg")
    velocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "velocity", "m/s")
    angularFrequency = carParameter{carVar}(0.0, "angular velocity", "rad/s")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Tire Force", "N")
    slipAngle = carParameter{carVar}(0.0, "Slip angle", "rad")
    slipRatio = carParameter{carVar}(0.0, "slip ratio", "-")
    maxSlipAngle = carParameter{carVar}(5/180*pi, "max slip angle", "rad")
    scalingForce = carParameter{carVar}(0.0, "max longitudinal force", "N")
    scalingForce.value = motor.torqueSpeedFunction(0.0) * gearbox.ratio.value / radius.value
    frictionCoefficient = carParameter{carVar}(1.0, "Friction Coefficient", "-", :sensitivity)
    rollingResistance = carParameter{carVar}(-0.010, "Rolling resistance coeff", "-", :sensitivity)
    brakingForce = carParameter{carVar}(0.0, "braking force", "N")
    # Pacejka Magic Formula coefficients
    B = 9.62   # stiffness factor
    C = 2.59  # shape factor
    E = 1  # curvature factor

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
        #static (non-dynamic) tire: angular frequency follows longitudinal velocity.
        #dynamic tire variants should override this and treat angularFrequency as a state.
        angularFrequency.value = velocityIn[1] / radius.value
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        α = slipAngle.value * 1
        D = forces.value[3] * frictionCoefficient.value

        if isnothing(optiModel)
            forces.value[2] = D * sin(C * atan(B * α - E * (B * α - atan(B * α))))
            forces.value[1] = inTorque/radius.value + brakingForce.value + forces.value[3] * rollingResistance.value
        else
            sf = scalingForce.value
            fy_expr = (D * sin(C * atan(B * α - E * (B * α - atan(B * α))))) / sf
            fx_expr = (inTorque/radius.value + brakingForce.value + forces.value[3] * rollingResistance.value) / sf
            fy_raw = @variable(optiModel)
            fx_raw = @variable(optiModel)
            @constraint(optiModel, fy_raw == fy_expr)
            @constraint(optiModel, fx_raw == fx_expr)
            forces.value[2] = fy_raw * sf
            forces.value[1] = fx_raw * sf
        end
    end

    function tireConstraints(model=nothing)
        μFz = frictionCoefficient.value * forces.value[3]
        fx = forces.value[1] / μFz
        fy = forces.value[2] / μFz
        lessContraint(fx^2 + fy^2, 1.1, model) #crucial to remove oscilatory behavior with bus, should look deeper into tire model suitable for optimisation
    end

    function setupObservables(ax::Axis)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        rect_obs = Observable(dummy)
        poly!(ax, rect_obs; color=(:black, 0.3), strokecolor=(:gray30, 0.5), strokewidth=1)
        ellipse_obs = Observable([Point2f(0, 0) for _ in _ELLIPSE_ANGLES])
        lines!(ax, ellipse_obs; color=(:gray, 0.5), linewidth=1, linestyle=:dash)
        force_pos_obs = Observable([Point2f(0, 0)])
        force_dir_obs = Observable([Point2f(0, 0)])
        arrows2d!(ax, force_pos_obs, force_dir_obs; color=:green, shaftwidth=4, tipwidth=8, tiplength=8)
        return (rect=rect_obs, ellipse=ellipse_obs, force_pos=force_pos_obs, force_dir=force_dir_obs)
    end

    function updateObservables(obs::NamedTuple, wx::Float64, wy::Float64, θ::Float64, Fz_static::Float64)
        r = radius.value
        s = 4r / Fz_static
        R = _rotmat2d(θ)
        obs.rect[] = _rect_points(wx, wy, 2r, width.value * 0.8, θ)
        Fmax_s = abs(forces.value[3]) * frictionCoefficient.value * s
        obs.ellipse[] = [Point2f(R * [Fmax_s * cos(a), Fmax_s * sin(a)] + [wx, wy]) for a in _ELLIPSE_ANGLES]
        f_global = R * [forces.value[1], forces.value[2]] .* s
        obs.force_pos[] = [Point2f(wx, wy)]
        obs.force_dir[] = [Point2f(f_global[1], f_global[2])]
    end

    tire = Tire(
        radius,
        width,
        inertia,
        mass,
        velocity,
        angularFrequency,
        forces,
        slipAngle,
        slipRatio,
        compute,
        tireConstraints,
        setVelocity,
        maxSlipAngle,
        scalingForce,
        frictionCoefficient,
        rollingResistance,
        brakingForce,
        setupObservables,
        updateObservables,
    )
end

function createR20_perantoni(motor::Motor, gearbox::Gearbox)
    # Perantoni & Limebeer (2014), "Optimal Control for a Formula One Car with Variable Parameters",
    # Appendix A + Table 3. Load-sensitive peak friction with linear interpolation between
    # two reference vertical loads, combined slip via normalised slip ρ and Pacejka-style shape factor.
    radius = carParameter{carVar}(0.205, "tire radius", "m")
    width = carParameter{carVar}(0.3, "tire width", "m")
    inertia = carParameter{carVar}(0.3, "tire inertia", "kg m^2")
    mass = carParameter{carVar}(1.0, "tire mass", "kg")
    velocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "velocity", "m/s")
    angularFrequency = carParameter{carVar}(0.0, "angular velocity", "rad/s")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Tire Force", "N")
    slipAngle = carParameter{carVar}(0.0, "Slip angle", "rad")
    slipRatio = carParameter{carVar}(0.0, "slip ratio", "-")
    maxSlipAngle = carParameter{carVar}(5/180*pi, "max slip angle (peak Fy at Fz1)", "rad")
    scalingForce = carParameter{carVar}(0.0, "max longitudinal force", "N")
    scalingForce.value = motor.torqueSpeedFunction(0.0) * gearbox.ratio.value / radius.value
    frictionCoefficient = carParameter{carVar}(1.0, "Global friction scaler", "-", :sensitivity)
    rollingResistance = carParameter{carVar}(0.010, "Rolling resistance coeff", "-", :sensitivity)
    brakingForce = carParameter{carVar}(0.0, "braking force", "N")

    # Perantoni-Limebeer structure scaled to FSAE R20 (Hoosier-class) operating range.
    # Reference loads match per-tyre static..peak transfer for ~250 kg car.
    # Peak coefficients chosen so pure-cornering Fy/Fz ≈ 1.0 (peak factor sin(Q·atan(S·1)) ≈ 0.965).
    Fz1, Fz2 = 500.0, 1500.0       # reference loads (N) — FSAE per-tyre range
    μx1, μx2 = 1.10, 0.95          # peak longitudinal friction at Fz1, Fz2
    μy1, μy2 = 1.10, 0.95          # peak lateral friction at Fz1, Fz2
    κ1, κ2   = 0.11, 0.10          # slip coefficient at peak Fx
    α1, α2   = 5/180*pi, 4/180*pi  # slip angle at peak Fy (FSAE tyres peak earlier than F1)
    Qx, Qy   = 1.9, 1.9            # shape factors (unchanged — dimensionless)
    Sx = pi / (2 * atan(Qx))
    Sy = pi / (2 * atan(Qy))
    ρ_eps = 1e-6                   # regulariser for ρ→0

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
        angularFrequency.value = velocityIn[1] / radius.value
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        Fz = forces.value[3]
        μx_max = (Fz - Fz1) * (μx2 - μx1) / (Fz2 - Fz1) + μx1   # eq. A.3
        μy_max = (Fz - Fz1) * (μy2 - μy1) / (Fz2 - Fz1) + μy1   # eq. A.4
        κ_peak = (Fz - Fz1) * (κ2  - κ1 ) / (Fz2 - Fz1) + κ1    # eq. A.5
        α_peak = (Fz - Fz1) * (α2  - α1 ) / (Fz2 - Fz1) + α1    # eq. A.6
        μx_max *= frictionCoefficient.value
        μy_max *= frictionCoefficient.value

        κ = slipRatio.value
        α = slipAngle.value
        κn = κ / κ_peak                                          # eq. A.7
        αn = α / α_peak                                          # eq. A.8
        ρ  = sqrt(αn^2 + κn^2 + ρ_eps)                           # eq. A.10
        μx = μx_max * sin(Qx * atan(Sx * ρ))                     # eq. A.11
        μy = μy_max * sin(Qy * atan(Sy * ρ))                     # eq. A.12
        Fy_tire = μy * Fz * αn / ρ                               # eq. A.16
        Fx_drive = inTorque/radius.value + brakingForce.value + Fz * rollingResistance.value

        if isnothing(optiModel)
            forces.value[2] = Fy_tire
            forces.value[1] = Fx_drive
        else
            sf = scalingForce.value
            fx_raw = @variable(optiModel, lower_bound=-2.0, upper_bound=2.0, start=0.0)
            fy_raw = @variable(optiModel, lower_bound=-2.0, upper_bound=2.0, start=0.0)
            @constraint(optiModel, fy_raw * sf == Fy_tire)
            @constraint(optiModel, fx_raw * sf == Fx_drive)
            forces.value[2] = fy_raw * sf
            forces.value[1] = fx_raw * sf
        end
    end

    function tireConstraints(model=nothing)
        # Load-sensitive friction ellipse using peak coefficients at current Fz.
        Fz = forces.value[3]
        μx_max = ((Fz - Fz1) * (μx2 - μx1) / (Fz2 - Fz1) + μx1) * frictionCoefficient.value
        μy_max = ((Fz - Fz1) * (μy2 - μy1) / (Fz2 - Fz1) + μy1) * frictionCoefficient.value
        fx = forces.value[1] / (μx_max * Fz)
        fy = forces.value[2] / (μy_max * Fz)
        lessContraint(fx^2 + fy^2, 1.0, model)
    end

    function setupObservables(ax::Axis)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        rect_obs = Observable(dummy)
        poly!(ax, rect_obs; color=(:black, 0.3), strokecolor=(:gray30, 0.5), strokewidth=1)
        ellipse_obs = Observable([Point2f(0, 0) for _ in _ELLIPSE_ANGLES])
        lines!(ax, ellipse_obs; color=(:gray, 0.5), linewidth=1, linestyle=:dash)
        force_pos_obs = Observable([Point2f(0, 0)])
        force_dir_obs = Observable([Point2f(0, 0)])
        arrows2d!(ax, force_pos_obs, force_dir_obs; color=:green, shaftwidth=4, tipwidth=8, tiplength=8)
        return (rect=rect_obs, ellipse=ellipse_obs, force_pos=force_pos_obs, force_dir=force_dir_obs)
    end

    function updateObservables(obs::NamedTuple, wx::Float64, wy::Float64, θ::Float64, Fz_static::Float64)
        r = radius.value
        s = 4r / Fz_static
        R = _rotmat2d(θ)
        obs.rect[] = _rect_points(wx, wy, 2r, width.value * 0.8, θ)
        Fz = forces.value[3]
        μ_avg = 0.5 * (((Fz - Fz1)*(μx2-μx1)/(Fz2-Fz1)+μx1) + ((Fz - Fz1)*(μy2-μy1)/(Fz2-Fz1)+μy1)) * frictionCoefficient.value
        Fmax_s = abs(Fz) * μ_avg * s
        obs.ellipse[] = [Point2f(R * [Fmax_s * cos(a), Fmax_s * sin(a)] + [wx, wy]) for a in _ELLIPSE_ANGLES]
        f_global = R * [forces.value[1], forces.value[2]] .* s
        obs.force_pos[] = [Point2f(wx, wy)]
        obs.force_dir[] = [Point2f(f_global[1], f_global[2])]
    end

    tire = Tire(
        radius,
        width,
        inertia,
        mass,
        velocity,
        angularFrequency,
        forces,
        slipAngle,
        slipRatio,
        compute,
        tireConstraints,
        setVelocity,
        maxSlipAngle,
        scalingForce,
        frictionCoefficient,
        rollingResistance,
        brakingForce,
        setupObservables,
        updateObservables,
    )
end

function create_FormulaE_pacejka(motor::Motor,gearbox::Gearbox)
    radius = carParameter{carVar}(0.33, "tire radius", "m")
    width = carParameter{carVar}(0.3, "tire width", "m")
    inertia = carParameter{carVar}(0.0, "tire inerti", "m")
    mass = carParameter{carVar}(0.0, "tire mass,wrong", "kg")
    velocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "velocity", "m/s")
    angularFrequency = carParameter{carVar}(0.0, "angular velocity", "rad/s")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Tire Force", "N")
    slipAngle = carParameter{carVar}(0.0, "Slip angle", "rad")
    slipRatio = carParameter{carVar}(0.0, "slip ratio", "-")
    maxSlipAngle = carParameter{carVar}(5/180*pi, "max slip angle", "rad")
    scalingForce = carParameter{carVar}(0.0, "max longitudinal force", "N")
    scalingForce.value = motor.torqueSpeedFunction(0.0) * gearbox.ratio.value / radius.value
    frictionCoefficient = carParameter{carVar}(1.0, "Friction Coefficient", "-", :sensitivity)
    rollingResistance = carParameter{carVar}(-0.010, "Rolling resistance coeff", "-")
    brakingForce = carParameter{carVar}(0.0, "braking force", "N")
    # Pacejka Magic Formula coefficients
    B = 9.62   # stiffness factor
    C = 2.59  # shape factor
    E = 1  # curvature factor

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
        #static (non-dynamic) tire: angular frequency follows longitudinal velocity.
        #dynamic tire variants should override this and treat angularFrequency as a state.
        angularFrequency.value = velocityIn[1] / radius.value
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        α = slipAngle.value * 1
        D = forces.value[3] * frictionCoefficient.value*(1 + -0.0813 * forces.value[3]/3000)*0.5
        forces.value[2] = D * sin(C * atan(B * α - E * (B * α - atan(B * α))))
        forces.value[1] = inTorque/radius.value + brakingForce.value + forces.value[3] * rollingResistance.value
    end

    function tireConstraints(model=nothing)
        μFz = frictionCoefficient.value * forces.value[3]
        fx = forces.value[1] / μFz
        fy = forces.value[2] / μFz
        lessContraint(fx^2 + fy^2, 1.0, model)
    end

    function setupObservables(ax::Axis)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        rect_obs = Observable(dummy)
        poly!(ax, rect_obs; color=(:black, 0.3), strokecolor=(:gray30, 0.5), strokewidth=1)
        ellipse_obs = Observable([Point2f(0, 0) for _ in _ELLIPSE_ANGLES])
        lines!(ax, ellipse_obs; color=(:gray, 0.5), linewidth=1, linestyle=:dash)
        force_pos_obs = Observable([Point2f(0, 0)])
        force_dir_obs = Observable([Point2f(0, 0)])
        arrows2d!(ax, force_pos_obs, force_dir_obs; color=:green, shaftwidth=4, tipwidth=8, tiplength=8)
        return (rect=rect_obs, ellipse=ellipse_obs, force_pos=force_pos_obs, force_dir=force_dir_obs)
    end

    function updateObservables(obs::NamedTuple, wx::Float64, wy::Float64, θ::Float64, Fz_static::Float64)
        r = radius.value
        s = 4r / Fz_static
        R = _rotmat2d(θ)
        obs.rect[] = _rect_points(wx, wy, 2r, width.value * 0.8, θ)
        Fmax_s = abs(forces.value[3]) * frictionCoefficient.value * s
        obs.ellipse[] = [Point2f(R * [Fmax_s * cos(a), Fmax_s * sin(a)] + [wx, wy]) for a in _ELLIPSE_ANGLES]
        f_global = R * [forces.value[1], forces.value[2]] .* s
        obs.force_pos[] = [Point2f(wx, wy)]
        obs.force_dir[] = [Point2f(f_global[1], f_global[2])]
    end

    tire = Tire(
        radius,
        width,
        inertia,
        mass,
        velocity,
        angularFrequency,
        forces,
        slipAngle,
        slipRatio,
        compute,
        tireConstraints,
        setVelocity,
        maxSlipAngle,
        scalingForce,
        frictionCoefficient,
        rollingResistance,
        brakingForce,
        setupObservables,
        updateObservables,
    )
end

