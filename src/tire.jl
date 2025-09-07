mutable struct  Tire
    radius::carParameter
    width::carParameter
    inertia::carParameter
    mass::carParameter
    velocity::carParameter
    angularFrequency::carParameter
    vertForce::carParameter
    longForce::carParameter
    latForce::carParameter
    slipAngle::carParameter
    slipRatio::carParameter
    tireFunction
    # where T
    #ire(radius, width, inertia, mass, angularFrequency, vertForce, longForce, latForce, slipAngle, slipRatio, tireFunction)
end
Base.show(io::IO, ::MIME"text/plain", obj::Tire) = prettyPrintComponent(io, obj)


#Pretty printing for tire
#function Base.show(io::IO, ::MIME"text/plain", obj::Tire)
#    T = typeof(obj)
#    println(io, "$(T)(")
#    # Calculate maximum field name length for alignment
#    max_width = maximum(length(string(field)) for field in fieldnames(T))
#    for field in fieldnames(T)
#        # Pad field names with spaces for alignment
#        field_str = rpad(string(field), max_width)
#        println(io, "  $(field_str) = ", getfield(obj, field))
#    end
#    print(io, ")")
#end




function createR20lin(maxTorque)
    radius = carParameter(0.205, "tire radius", "m")
    width = carParameter(0.3, "tire width, wrong", "m")
    inertia = carParameter(0.3, "tire width, wrong", "m")
    mass = carParameter(1.0, "tire mass,wrong", "kg")
    velocity = carParameter([0.0, 0.0, 0.0], "velocity", "m/s");
    angularFrequency = carParameter(0.0, "angular velocity", "rad/s")
    vertForce = carParameter(1000.0, "Vertical force on tire", "N")#palceholder value
    longForce = carParameter(0.0, "Longitudinal force from tire", "N")
    latForce = carParameter(0.0, "Lateral force from tire", "N")
    slipAngle = carParameter(0.0, "slip angle", "rad")
    slipRatio = carParameter(0.0, "slip ratio", "-")

    maxForce = maxTorque/radius.value

    function tireFunction(inTorque,optiModel=nothing)
        slipAngle.value = atan(velocity.value[1], velocity.value[2])
        latForce.value = slipAngle.value * vertForce.value
        longForce.value = inTorque
        if isnothing(optiModel)
            

        else
            @constraint(optiModel, (tire.latForce/maxForce)^2 + (tire.longForce/maxForce)^2 <= (tire.vertForce/maxForce)^2)

        end
        #tire.slipRatio = tire.angularFrequency * tire.radius / velocity[1]
        #tire.longForce = tire.slipRatio * vertForce
    end
    tire = Tire(
        radius,
        width,
        inertia,
        mass,
        velocity,
        angularFrequency,
        vertForce,
        longForce,
        latForce,
        slipAngle,
        slipRatio,
        tireFunction
    )
end