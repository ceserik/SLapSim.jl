mutable struct carParameter{T}
    value::T
    name::String
    unit::String
    size::Tuple{Int,Int}  # To store dimensions of the value (1,1) for scalar, (n,1) for vector, (n,m) for matrix
    role::Symbol  # :control, :state, :static parameter, :tunable parameter
end
function carParameter{T}(value::T, name::String, unit::String, role::Symbol=:static) where T
    if isa(value, AbstractArray)
        size_tuple = (length(value), 1)  # For vectors, treat as (n,1)
        if ndims(value) == 2
            size_tuple = size(value)  # For matrices, use actual dimensions
        end
    else
        size_tuple = (1, 1)  # For scalars
    end

    return carParameter{T}(value, name, unit, size_tuple, role)
end

# Convenience constructor for vector carParameters
carParameter{Vector{T}}(v::Vector{<:Real}, name::String, unit::String) where T = 
    carParameter{Vector{T}}(Vector{T}(v), name, unit)

# Pretty printing for carParameter
function Base.show(io::IO, p::carParameter)
    val = p.value
    if isa(val, AbstractArray)
        print(io, "$(p.name) = $(val) [$(p.unit)] ($(p.role))")
    else
        print(io, "$(p.name) = $(val) [$(p.unit)] ($(p.role))")
    end
end

function Base.show(io::IO, ::MIME"text/plain", p::carParameter)
    println(io, "carParameter: $(p.name)")
    println(io, "  Value: $(p.value)")
    println(io, "  Unit:  $(p.unit)")
    println(io, "  Size:  $(p.size)")
    print(io,   "  Role:  $(p.role)")
end

mutable struct CarParameters
    mass::carParameter{carVar}
    inertia::carParameter{carVar}
    motorForce::carParameter{carVar}
    CL::carParameter{carVar}
    CD::carParameter{carVar}
    velocity::carParameter{Vector{carVar}}
    angularVelocity::carParameter{Vector{carVar}}
    psi::carParameter{carVar}
    n::carParameter{carVar}
    powerLimit::carParameter{carVar}
    lateralForce::carParameter{carVar}
    lateralTransfer::carParameter{carVar}
    brakeBias::carParameter{carVar}
    nControls::carParameter{carVar}
    nStates::carParameter{carVar}
    s::carParameter{carVar}
end

# Pretty printing for CarParameters
function Base.show(io::IO, ::MIME"text/plain", car::CarParameters)
    println(io, "CarParameters:")
    for fname in fieldnames(CarParameters)
        p = getfield(car, fname)
        println(io, "  $(p.name) = $(p.value) [$(p.unit)] ($(p.role))")
    end
end