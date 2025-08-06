# Define a Parameter struct
mutable struct carParameter{T}
    value::Any  
    name::String
    unit::String
    
    # Constructor for numeric values
    carParameter(value::T, name::String, unit::String) where T = new{T}(value, name, unit)
end

mutable struct carParameters{T}
    mass::carParameter
    motorForce::carParameter
    CL::carParameter
    CD::carParameter
    vx::carParameter
    vy::carParameter
    psi::carParameter
    n::carParameter
    powerLimit::carParameter
    
end

# Constructor
function carParameters(
    mass::carParameter{T},
    motorForce::carParameter{T},
    CL::carParameter{T},
    CD::carParameter{T},
    vx::carParameter{T},
    vy::carParameter{T},
    psi::carParameter{T},
    n::carParameter{T},
    powerLimit::carParameter{T}) where T
    carParameters{T}(mass, motorForce, CL, CD, vx,vy,psi,n,powerLimit)
end

# Define the Car struct with parameters
mutable struct Car
    carFunction
    carParameters
    controlMapping
    stateMapping
    mapping

    # Inner constructor
    function Car(carFunction, carParameters, controlMapping = nothing, stateMapping = nothing, mapping = nothing)
        new(carFunction, carParameters, controlMapping,stateMapping,mapping)
    end
end


