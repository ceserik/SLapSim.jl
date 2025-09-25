mutable struct carParameter{T}
    value::T
    name::String
    unit::String
    size::Tuple{Int,Int}  # To store dimensions of the value (1,1) for scalar, (n,1) for vector, (n,m) for matrix
    
end
function carParameter{T}(value::T, name::String, unit::String) where T
    if isa(value, AbstractArray)
        size_tuple = (length(value), 1)  # For vectors, treat as (n,1)
        if ndims(value) == 2
            size_tuple = size(value)  # For matrices, use actual dimensions
        end
    else
        size_tuple = (1, 1)  # For scalars
    end
    
    return carParameter{T}(value, name, unit, size_tuple)
end

# Convenience constructor for vector carParameters
carParameter{Vector{T}}(v::Vector{<:Real}, name::String, unit::String) where T = 
    carParameter{Vector{T}}(Vector{T}(v), name, unit)

mutable struct CarParameters
    mass::carParameter{T} where T
    inertia::carParameter{T} where T
    motorForce::carParameter{T} where T
    CL::carParameter{T} where T
    CD::carParameter{T} where T
    velocity::carParameter{T} where T
    angularVelocity::carParameter{T} where T
    psi::carParameter{T} where T
    n::carParameter{T} where T
    powerLimit::carParameter{T} where T
    lateralForce::carParameter{T} where T
    nControls::carParameter{T} where T
    nStates::carParameter{T} where T
    s::carParameter{T} where T
end