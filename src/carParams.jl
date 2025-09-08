
mutable struct carParameter
    value::Any
    name::String
    unit::String
    size::Tuple{Int,Int}  # To store dimensions of the value (1,1) for scalar, (n,1) for vector, (n,m) for matrix
    
    # Constructor for scalar values
    function carParameter(value::Number, name::String, unit::String)
        new(value, name, unit, (1,1))
    end
    
    # Constructor for vector values
    function carParameter(value::AbstractVector, name::String, unit::String)
        new(value, name, unit, (length(value),1))
    end
    
    # Constructor for matrix values
    function carParameter(value::AbstractMatrix, name::String, unit::String)
        new(value, name, unit, size(value))
    end
end



mutable struct Drivetrain
    motors
    gearboxes
    tires
    accumulators
end

# Define the Car struct with parameters
mutable struct Car
    carFunction
    carParameters
    controlMapping
    stateMapping
    mapping

    # Inner constructor
    function Car(carFunction, carParameters, controlMapping=nothing, stateMapping=nothing, mapping=nothing)
        new(carFunction, carParameters, controlMapping, stateMapping, mapping)
    end
end

mutable struct Chassis
    mass::carParameter
end



function createCTU25chassis()
    mass = carParameter(180,"mass","kg")
    chassis = Chassis(
        mass
    )
    return chassis

end

mutable struct carParameters
    mass::carParameter
    motorForce::carParameter
    CL::carParameter
    CD::carParameter
    velocity::carParameter
    angularVelocity::carParameter
    psi::carParameter
    n::carParameter
    powerLimit::carParameter
    lateralForce::carParameter
    nControls::carParameter
end


function createCTU25_2D()
    mass = carParameter(280.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    lateralForce = carParameter(0.0, "lateral Force", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    velocity = carParameter([0,0,0],"Speed X","m/s")
    powerLimit = carParameter(80000.0,"PowerLimit","W")
    psi = carParameter(0.0,"heading","rad")
    n = carParameter(0.0,"Distance from centerline","m")
    nControls = carParameter(2.0,"number of controlled parameters","-")
    p = carParameters(
        mass,
        motorForce,
        CL,
        CD,
        velocity,
        angularVelocity,
        psi,
        n,
        powerLimit,
        lateralForce,
        nControls
    )

    function controlMapping(input, controls)
        input.motorForce.value = controls[1]  
        #paramcopy.CL.value = controls[2]
        return input
    end

    function mapping(car,u,x)
        car.carParameters.motorForce.value = u[1]
        car.carParameters.lateralForce.value = u[2]

        car.carParameters.vx.value = x[1]
        car.carParameters.vy.value = x[2]

    end

    function stateMapping(input,states)
        input.vx.value = states[1]
        #map psi to track heading here, maybe input should be also track mapping should be unified
        return input
    end
    car = Car(
        massPointCar,
        p,
        controlMapping,
        stateMapping,
        mapping
    )

    return car
end

mutable struct Car2
    drivetrain
    aero
    suspension
    chassis
    wheelAssemblies
    carParameters
end

    
function createSimplestSingleTrack()

    gearboxFront = createCTU25gearbox()
    gearboxRear = createCTU25gearbox()

    motorFront = createFischerMotor()
    motorRear = createFischerMotor()

    #get max motor torque for scaling
    maxMotorTorqueFront = motorFront.torqueSpeedFunction(0) * gearboxFront.ratio.value
    maxMotorTorqueRear = motorFront.torqueSpeedFunction(0) * gearboxRear.ratio.value

    tireFront = createR20lin(maxMotorTorqueFront)
    tireRear = createR20lin(maxMotorTorqueRear)


    drivetrain = Drivetrain(
        [motorFront,motorRear],
        [gearboxFront, gearboxRear],
        [tireFront,tireRear],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createDummySuspension()
    wheelAssemblyFront = createBasicWheelAssembly([1.520/2, 0, 0]) # wheelbase musi byt parameter!!!!
    wheelAssemblyRear = createBasicWheelAssembly([-1.520/2, 0, 0])
    chassis = createCTU25chassis()

    velocity = carParameter([0.0, 0.0, 0.0], "translational velocity", "m/s");
    angularVelocity = carParameter([0.0, 0.0, 0.0], "angular velocity", "rad/s");
    
    mass = carParameter(280.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    lateralForce = carParameter(0.0, "lateral Force", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    velocity = carParameter([0,0,0],"Speed X","m/s")
    powerLimit = carParameter(80000.0,"PowerLimit","W")
    psi = carParameter(0.0,"heading","rad")
    n = carParameter(0.0,"Distance from centerline","m")
    nControls = carParameter(2.0,"number of controlled parameters","-")
    
    p = carParameters(
        mass,
        motorForce,
        CL,
        CD,
        velocity,
        angularVelocity,
        psi,
        n,
        powerLimit,
        lateralForce,
        nControls
    )

    afto = Car2(
        drivetrain,
        aero,
        suspension,
        chassis,
        [wheelAssemblyFront,wheelAssemblyRear],
        p
    )
    return afto
end


"""
Generic pretty printing for car components.
Usage: import this function and add:
    Base.show(io::IO, ::MIME"text/plain", obj::YourType) = prettyPrintComponent(io, obj)
"""
function prettyPrintComponent(io::IO, obj)
    # List of types that should get pretty printing
    pretty_print_types = [Tire, Motor, Gearbox, Chassis, Drivetrain, Car2,WheelAssembly]
    
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

