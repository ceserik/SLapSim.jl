



#
#
#function createCTU25_2D()
#    mass = carParameter(280.0, "Mass", "kg")
#    inertia = carParameter(100.0, "Inertia", "kg*m^2")
#    motorForce = carParameter(1000.0, "motorForce", "N")
#    lateralForce = carParameter(0.0, "lateral Force", "N")
#    CL = carParameter(5.0, "Lift Coefficient", "-")
#    CD = carParameter(2.0, "Drag Coefficient", "-")
#    velocity = carParameter([0,0,0],"Speed X","m/s")
#    powerLimit = carParameter(80000.0,"PowerLimit","W")
#    psi = carParameter(0.0,"heading","rad")
#    n = carParameter(0.0,"Distance from centerline","m")
#    nControls = carParameter(2.0,"number of controlled parameters","-")
#    p = carParameters(
#        mass,
#        inertia,
#        motorForce,
#        CL,
#        CD,
#        velocity,
#        angularVelocity,
#        psi,
#        n,
#        powerLimit,
#        lateralForce,
#        nControls
#    )
#
#    function controlMapping(input, controls)
#        input.motorForce.value = controls[1]  
#        #paramcopy.CL.value = controls[2]
#        return input
#    end
#
#    function mapping(car,u,x)
#        car.carParameters.motorForce.value = u[1]
#        car.carParameters.lateralForce.value = u[2]
#
#        car.carParameters.vx.value = x[1]
#        car.carParameters.vy.value = x[2]
#
#    end
#
#    function stateMapping(input,states)
#        input.vx.value = states[1]
#        #map psi to track heading here, maybe input should be also track mapping should be unified
#        return input
#    end
#    car = Car(
#        massPointCar,
#        p,
#        controlMapping,
#        stateMapping,
#        mapping
#    )
#
#    return car
#end
#

    


"""
Generic pretty printing for car components.
Usage: import this function and add:
    Base.show(io::IO, ::MIME"text/plain", obj::YourType) = prettyPrintComponent(io, obj)
"""
function prettyPrintComponent(io::IO, obj)
    # List of types that should get pretty printing
    pretty_print_types = [Tire, Motor, Gearbox, Chassis, Drivetrain, Car2, WheelAssembly]
    
    # Handle types that shouldn't get pretty printing
    if !any(T -> obj isa T, pretty_print_types)
        print(io, obj)
        return
    end
    
    # Handle car components with pretty printing
    T = typeof(obj)
    println(io, "$(T)(")
    
    if fieldcount(T) > 0
        max_width = maximum(length(string(field)) for field in fieldnames(T))
        for field in fieldnames(T)
            field_str = rpad(string(field), max_width)
            println(io, "  $(field_str) = ", getfield(obj, field))
        end
    end
    
    print(io, ")")
end


print("carParams.jl great success\n")

