mutable struct carParameter{T}
    value::T
    name::String
    unit::String
    size::Tuple{Int,Int}  # To store dimensions of the value (1,1) for scalar, (n,1) for vector, (n,m) for matrix
    role::Symbol  # :control, :state, :static parameter, :tunable parameter
    limits::Vector{Float64}
end
function carParameter{T}(value::T, name::String, unit::String, role::Symbol=:static,limits::Vector{Float64} = [9999.0,-9999.0]) where T
    if isa(value, AbstractArray)
        size_tuple = (length(value), 1)  # For vectors, treat as (n,1)
        if ndims(value) == 2
            size_tuple = size(value)  # For matrices, use actual dimensions
        end
    else
        size_tuple = (1, 1)  # For scalars
    end

    return carParameter{T}(value, name, unit, size_tuple, role,limits)
end

# Convenience constructor for vector carParameters
carParameter{Vector{T}}(v::Vector{<:Real}, name::String, unit::String, args...) where T =
    carParameter{Vector{T}}(Vector{T}(v), name, unit, args...)

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
    println(io, "  Value:  $(p.value)")
    println(io, "  Unit:   $(p.unit)")
    println(io, "  Size:   $(p.size)")
    println(io, "  Role:   $(p.role)")
    print(io,   "  Limits: $(p.limits)")
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

# --- Variable descriptors for optimizer mapping ---

struct VarEntry
    name::String
    targets::Vector{Pair}  # carParameter => vector_index (0 = scalar)
    lb::Float64
    ub::Float64
end

# Auto-read limits from the first target's carParameter
function VarEntry(name::String, targets::Vector{<:Pair})
    param = first(targets).first
    VarEntry(name, collect(Pair, targets), param.limits[1], param.limits[2])
end

# Explicit bounds with typed pair vector
function VarEntry(name::String, targets::Vector{<:Pair}, lb::Real, ub::Real)
    VarEntry(name, collect(Pair, targets), Float64(lb), Float64(ub))
end

function apply_mapping!(desc::Vector{VarEntry}, values::AbstractVector)
    for (i, entry) in enumerate(desc)
        for (param, idx) in entry.targets
            if idx == 0
                param.value = values[i]
            else
                param.value[idx] = values[i]
            end
        end
    end
end

get_bounds(desc::Vector{VarEntry}) = ([e.lb for e in desc], [e.ub for e in desc])

function get_scales(desc::Vector{VarEntry})
    return [max(abs(e.lb), abs(e.ub), 1.0) for e in desc]
end