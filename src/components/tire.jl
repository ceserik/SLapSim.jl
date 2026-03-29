using Revise
mutable struct  Tire{F1,F2,F3}
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
end
Base.show(io::IO, ::MIME"text/plain", obj::Tire) = prettyPrintComponent(io, obj)

function createR20lin_double(motor,gearbox)
    radius = carParameter{carVar}(0.205, "tire radius", "m")
    width = carParameter{carVar}(0.3, "tire width, wrong", "m")
    inertia = carParameter{carVar}(0.3, "tire width, wrong", "m")
    mass = carParameter{carVar}(1.0, "tire mass,wrong", "kg")
    velocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "velocity", "m/s")
    angularFrequency = carParameter{carVar}(0.0, "angular velocity", "rad/s")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Forces from tire x y z", "N")
    slipAngle = carParameter{carVar}(0.0, "slip angle", "rad")
    slipRatio = carParameter{carVar}(0.0, "slip ratio", "-")
    maxSlipAngle = carParameter{carVar}(5/180*pi, "max slip angle", "rad")
    scalingForce = carParameter{carVar}(0.0, "max longitudinal force", "N")
    scalingForce.value = motor.torqueSpeedFunction(0.0) * gearbox.ratio.value / radius.value

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        forces.value[2] = 2*slipAngle.value/maxSlipAngle.value * forces.value[3]
        forces.value[1] = inTorque/radius.value
    end

    function tireConstraints(model=nothing)
        lessContraint((forces.value[2]/scalingForce.value)^2 + (forces.value[1]/scalingForce.value)^2, (forces.value[3]/scalingForce.value)^2, model)
        lessContraint(slipAngle.value/maxSlipAngle.value, 5/180*pi/maxSlipAngle.value, model)
        greaterContraint(slipAngle.value/maxSlipAngle.value, -5/180*pi/maxSlipAngle.value, model)
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
    )
    return tire
end

function createR20lin(motor,gearbox)
    radius = carParameter{carVar}(0.205, "tire radius", "m")
    width = carParameter{carVar}(0.3, "tire width, wrong", "m")
    inertia = carParameter{carVar}(0.3, "tire width, wrong", "m")
    mass = carParameter{carVar}(1.0, "tire mass,wrong", "kg")
    velocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "velocity", "m/s")
    angularFrequency = carParameter{carVar}(0.0, "angular velocity", "rad/s")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Forces from tire x y z", "N")
    slipAngle = carParameter{carVar}(0.0, "slip angle", "rad")
    slipRatio = carParameter{carVar}(0.0, "slip ratio", "-")
    maxSlipAngle = carParameter{carVar}(5/180*pi, "max slip angle", "rad")
    scalingForce = carParameter{carVar}(0.0, "max longitudinal force", "N")
    scalingForce.value = motor.torqueSpeedFunction(0.0) * gearbox.ratio.value / radius.value

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        forces.value[2] = 1*slipAngle.value/maxSlipAngle.value * forces.value[3]
        forces.value[1] = inTorque/radius.value
    end

    function tireConstraints(model=nothing)
        lessContraint((forces.value[2]/scalingForce.value)^2 + (forces.value[1]/scalingForce.value)^2, (forces.value[3]/scalingForce.value)^2, model)
        lessContraint(slipAngle.value/maxSlipAngle.value, 5/180*pi/maxSlipAngle.value, model)
        greaterContraint(slipAngle.value/maxSlipAngle.value, -5/180*pi/maxSlipAngle.value, model)
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
    )
end

"""
    draw!(ax, tire::Tire, wx, wy, θ; force_scale=1/1000)

Draw a tire as a filled rectangle at world position (wx,wy) with total orientation θ
(car heading + steering angle). Also draws the friction ellipse and resulting force vector.
`force_scale` converts force [N] to drawing units [m].
"""
function draw!(ax, tire::Tire, wx, wy, θ; force_scale=1/1000)
    r = tire.radius.value
    w = tire.width.value * 0.8
    c = _rect_corners(wx, wy, 2r, w, θ)
    poly!(ax, Point2f.(eachrow(c)); color=:black, strokecolor=:gray30, strokewidth=1)

    R = _rotmat2d(θ)
    Fz = tire.forces.value[3]

    # Friction ellipse (circle with radius Fz, scaled)
    n_pts = 40
    ellipse_r = abs(Fz) * force_scale
    angles = range(0, 2π, length=n_pts + 1)
    ellipse_pts = [Point2f(R * [ellipse_r * cos(a); ellipse_r * sin(a)] .+ [wx, wy]) for a in angles]
    lines!(ax, ellipse_pts; color=(:gray, 0.5), linewidth=1, linestyle=:dash)

    # Force vector (Fx, Fy in tire frame → rotated to global)
    Fx = tire.forces.value[1]
    Fy = tire.forces.value[2]
    f_global = R * [Fx; Fy] .* force_scale
    if norm(f_global) > 1e-6
        arrows2d!(ax, [wx], [wy], [f_global[1]], [f_global[2]];
                  color=:green, shaftwidth=1.5, tipwidth=6, tiplength=6)
    end
end

const _N_ELLIPSE_PTS = 41
const _ELLIPSE_ANGLES = range(0, 2π, length=_N_ELLIPSE_PTS)

function setup_observables!(ax, tire::Tire)
    dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
    rect_obs = Observable(dummy)
    poly!(ax, rect_obs; color=:black, strokecolor=:gray30, strokewidth=1)

    ellipse_obs = Observable([Point2f(0, 0) for _ in _ELLIPSE_ANGLES])
    lines!(ax, ellipse_obs; color=(:gray, 0.5), linewidth=1, linestyle=:dash)

    force_pos_obs = Observable([Point2f(0, 0)])
    force_dir_obs = Observable([Point2f(0, 0)])
    arrows2d!(ax, force_pos_obs, force_dir_obs; color=:green, shaftwidth=1.5, tipwidth=6, tiplength=6)

    return (rect=rect_obs, ellipse=ellipse_obs, force_pos=force_pos_obs, force_dir=force_dir_obs)
end

function update_observables!(obs, tire::Tire, wx, wy, θ; force_scale=1/1000)
    r = tire.radius.value
    w = tire.width.value * 0.8
    obs.rect[] = _rect_points(wx, wy, 2r, w, θ)

    R = _rotmat2d(θ)
    Fz = tire.forces.value[3]
    ellipse_r = abs(Fz) * force_scale
    obs.ellipse[] = [Point2f(R * [ellipse_r * cos(a); ellipse_r * sin(a)] .+ [wx, wy]) for a in _ELLIPSE_ANGLES]

    Fx = tire.forces.value[1]
    Fy = tire.forces.value[2]
    f_global = R * [Fx; Fy] .* force_scale
    obs.force_pos[] = [Point2f(wx, wy)]
    obs.force_dir[] = [Point2f(f_global[1], f_global[2])]
end