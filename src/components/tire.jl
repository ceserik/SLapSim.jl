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
    draw!(ax, tire::Tire, wx, wy, θ)

Draw a tire as a filled rectangle at world position (wx,wy) with total orientation θ
(car heading + steering angle).
"""
function draw!(ax, tire::Tire, wx, wy, θ)
    r = tire.radius.value
    w = tire.width.value * 0.8
    c = _rect_corners(wx, wy, 2r, w, θ)
    poly!(ax, Point2f.(eachrow(c)); color=:black, strokecolor=:gray30, strokewidth=1)
end