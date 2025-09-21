mutable struct  Tire
    radius::carParameter{carVar}
    width::carParameter{carVar}
    inertia::carParameter{carVar}
    mass::carParameter{carVar}
    velocity::carParameter{Vector{carVar}}
    angularFrequency::carParameter{carVar}
    forces::carParameter{Vector{carVar}}
    slipAngle::carParameter{carVar}
    slipRatio::carParameter{carVar}
    tireFunction::Function
    
end
Base.show(io::IO, ::MIME"text/plain", obj::Tire) = prettyPrintComponent(io, obj)



function createR20lin(maxTorque::Float64)
    radius = carParameter{carVar}(0.205, "tire radius", "m")
    width = carParameter{carVar}(0.3, "tire width, wrong", "m")
    inertia = carParameter{carVar}(0.3, "tire width, wrong", "m")
    mass = carParameter{carVar}(1.0, "tire mass,wrong", "kg")
    velocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "velocity", "m/s");
    angularFrequency = carParameter{carVar}(0.0, "angular velocity", "rad/s")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Forces from tire x y z", "N")
    slipAngle = carParameter{carVar}(0.0, "slip angle", "rad")
    slipRatio = carParameter{carVar}(0.0, "slip ratio", "-")

    maxForce = maxTorque/radius.value

    function tireFunction(inTorque::carVar,optiModel::JuMP.Model=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        forces.value[2] = slipAngle.value * forces.value[3]
        forces.value[1] = inTorque/radius.value
        if isnothing(optiModel)

        else
            @constraint(optiModel, (tire.forces.value[2]/maxForce)^2 + (tire.forces.value[1]/maxForce)^2 <= (tire.forces.value[3]/maxForce)^2)
            @constraint(optiModel, slipAngle.value <=  5*0.08726646259971647)
            @constraint(optiModel, slipAngle.value >= -5*0.08726646259971647)

        end
        #tire.slipRatio = tire.angularFrequency * tire.radius / velocity[1]
        #tire.longForce = tire.slipRatio * 
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
        tireFunction
    )
end