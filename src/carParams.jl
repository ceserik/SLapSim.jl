
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


mutable struct carParameters
    mass::carParameter
    motorForce::carParameter
    CL::carParameter
    CD::carParameter
    vx::carParameter
    vy::carParameter
    psi::carParameter
    n::carParameter
    powerLimit::carParameter
    lateralForce::carParameter
    nControls::carParameter
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


mutable struct Car2
    velocity
    angularVelocity
    drivetrain
    aero
    suspension
    chassis
    wheelAssemblies

end

function createCTU25chassis()
    mass = carParameter(180,"mass","kg")
    chassis = Chassis(
        mass
    )
    return chassis

end

function createCTU25_2D()
    mass = carParameter(280.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    lateralForce = carParameter(0.0, "lateral Force", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    vx = carParameter(15.0,"Speed X","m/s")
    powerLimit = carParameter(80000.0,"PowerLimit","W")
    vy = carParameter(0.0,"Speed Y","m/s")
    psi = carParameter(0.0,"heading","rad")
    n = carParameter(0.0,"Distance from centerline","m")
    nControls = carParameter(2.0,"number of controlled parameters","-")
    p = carParameters(
        mass,
        motorForce,
        CL,
        CD,
        vx,
        vy,
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
    


function createSimplestSingleTrack()
    tireFront = createR20lin(500)
    tireRear = createR20lin(500)
    
    gearboxFront = createCTU25gearbox()
    gearboxRear = createCTU25gearbox()

    motorFront = createFischerMotor()
    motorRear = createFischerMotor()

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

    afto = Car2(
        velocity,
        angularVelocity,
        drivetrain,
        aero,
        suspension,
        chassis,
        [wheelAssemblyFront,wheelAssemblyRear]
    )
    return afto
end
print("carParams.jl great success\n")
