using Revise
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
    setupObservables::F4
    updateObservables::F5
end
Base.show(io::IO, ::MIME"text/plain", obj::Tire) = prettyPrintComponent(io, obj)

function createR20lin(motor,gearbox)
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

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        forces.value[2] = 1*slipAngle.value/maxSlipAngle.value * forces.value[3] * frictionCoefficient.value
        forces.value[1] = inTorque/radius.value
    end

    function tireConstraints(model=nothing)
        lessContraint((forces.value[2]/scalingForce.value)^2 + (forces.value[1]/scalingForce.value)^2, (frictionCoefficient.value * forces.value[3]/scalingForce.value)^2, model)
        #lessContraint(slipAngle.value/maxSlipAngle.value, maxSlipAngle.value/maxSlipAngle.value, model)
        #greaterContraint(slipAngle.value/maxSlipAngle.value, -maxSlipAngle.value/maxSlipAngle.value, model)
    end

    function setupObservables(ax)
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

    function updateObservables(obs, wx, wy, θ, Fz_static)
        r = radius.value
        s = 4r / Fz_static
        R = _rotmat2d(θ)
        obs.rect[] = _rect_points(wx, wy, 2r, width.value * 0.8, θ)
        obs.ellipse[] = [Point2f(R * [abs(forces.value[3]) * s * cos(a); abs(forces.value[3]) * s * sin(a)] .+ [wx, wy]) for a in _ELLIPSE_ANGLES]
        f_global = R * forces.value[1:2] .* s
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
        setupObservables,
        updateObservables,
    )
end

const _N_ELLIPSE_PTS = 41
const _ELLIPSE_ANGLES = range(0, 2π, length=_N_ELLIPSE_PTS)

