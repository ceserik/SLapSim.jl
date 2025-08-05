# Define a Parameter struct
mutable struct carParameter{T}
    value::Any  # Changed back to Any to allow both numeric values and optimization variables
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
end

# Constructor
function carParameters(mass::carParameter{T}, motorForce::carParameter{T}, CL::carParameter{T}, CD::carParameter{T}) where T
    carParameters{T}(mass, motorForce, CL, CD)
end

# Define the Car struct with parameters
mutable struct Car
    carFunction
    carParameters
    controlMapping

    # Inner constructor
    function Car(carFunction, carParameters, controlMapping = nothing)
        new(carFunction, carParameters, controlMapping)
    end
end


