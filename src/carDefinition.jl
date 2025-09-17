
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
    inertia::carParameter
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

mutable struct Car2
    drivetrain
    aero
    suspension
    chassis
    wheelAssemblies
    carParameters
    carFunction
    controlMapping::Function
    stateMapping::Function
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

