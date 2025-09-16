mutable struct  Tire
    radius::carParameter
    width::carParameter
    inertia::carParameter
    mass::carParameter
    velocity::carParameter
    angularFrequency::carParameter
    forces::carParameter
    slipAngle::carParameter
    slipRatio::carParameter
    tireFunction
    
end
Base.show(io::IO, ::MIME"text/plain", obj::Tire) = prettyPrintComponent(io, obj)



function createR20lin(maxTorque)
    radius = carParameter(0.205, "tire radius", "m")
    width = carParameter(0.3, "tire width, wrong", "m")
    inertia = carParameter(0.3, "tire width, wrong", "m")
    mass = carParameter(1.0, "tire mass,wrong", "kg")
    velocity = carParameter([0.0, 0.0, 0.0], "velocity", "m/s");
    angularFrequency = carParameter(0.0, "angular velocity", "rad/s")
    forces = carParameter([0.0, 0.0, 0.0], "Forces from tire x y z", "N")
    slipAngle = carParameter(0.0, "slip angle", "rad")
    slipRatio = carParameter(0.0, "slip ratio", "-")

    maxForce = maxTorque/radius.value

    function tireFunction(inTorque,optiModel=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        forces.value[2] = slipAngle.value * forces.value[3]
        forces.value[1] = inTorque/radius.value
        if isnothing(optiModel)

        else
            @constraint(optiModel, (tire.latForce/maxForce)^2 + (tire.longForce/maxForce)^2 <= (tire.force[3].value/maxForce)^2)
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