
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
    frictionCoefficient = carParameter{carVar}(1.0, "Friction Coefficient", "-", :tunable)
    rollingResistance = carParameter{carVar}(0.010, "Rolling resistance coeff", "-", :tunable)
    brakingForce = carParameter{carVar}(0.0, "braking force", "N")

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        forces.value[2] = 1*slipAngle.value/maxSlipAngle.value * forces.value[3] * frictionCoefficient.value
        forces.value[1] = inTorque/radius.value
    end

    function tireConstraints(model=nothing)
        Fx_s = forces.value[1]/scalingForce.value
        Fy_s = forces.value[2]/scalingForce.value
        Fz_s = frictionCoefficient.value * forces.value[3]/scalingForce.value
        if isnothing(model)
            lessContraint(Fx_s^2 + Fy_s^2, Fz_s^2, nothing)
        else
            ε = @variable(model, lower_bound=0.0)
            @constraint(model, Fx_s^2 + Fy_s^2 + ε == Fz_s^2)
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
    frictionCoefficient = carParameter{carVar}(1.0, "Friction Coefficient", "-", :tunable)
    rollingResistance = carParameter{carVar}(0.010, "Rolling resistance coeff", "-", :tunable)
    brakingForce = carParameter{carVar}(0.0, "braking force", "N")
    # Pacejka Magic Formula coefficients
    B = 5   # stiffness factor
    C = 2.59  # shape factor
    E = 1  # curvature factor

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        α = slipAngle.value * 1
        D = forces.value[3] * frictionCoefficient.value
        forces.value[2] = D * sin(C * atan(B * α - E * (B * α - atan(B * α))))
        forces.value[1] = inTorque/radius.value
    end

    function tireConstraints(model=nothing)
        Fx_s = forces.value[1]/scalingForce.value
        Fy_s = forces.value[2]/scalingForce.value
        Fz_s = frictionCoefficient.value * forces.value[3]/scalingForce.value
        if isnothing(model)
            lessContraint(Fx_s^2 + Fy_s^2, Fz_s^2, nothing)
        else
            ε = @variable(model, lower_bound=0.0)
            @constraint(model, Fx_s^2 + Fy_s^2 + ε == Fz_s^2)
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



