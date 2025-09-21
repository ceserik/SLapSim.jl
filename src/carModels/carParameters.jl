mutable struct carParameter{T}
    value::T
    name::String
    unit::String
    size::Tuple{Int,Int}  # To store dimensions of the value (1,1) for scalar, (n,1) for vector, (n,m) for matrix
    
#    # Constructor for scalar values
#    function carParameter(value::Number, name::String, unit::String)
#        new(value, name, unit, (1,1))
#    end
#    
#    # Constructor for vector values
#    function carParameter(value::AbstractVector, name::String, unit::String)
#        new(value, name, unit, (length(value),1))
#    end
#    
#    # Constructor for matrix values
#    function carParameter(value::AbstractMatrix, name::String, unit::String)
#        new(value, name, unit, size(value))
#    end
end


# outer constructors that produce concrete parameter types
carParameter(value::Number, name::String, unit::String) = carParameter{Float64}(float(value), name, unit, (1,1))
carParameter(value::AbstractVector{<:Number}, name::String, unit::String) = carParameter{Vector{Float64}}(Float64.(collect(value)), name, unit, (length(value),1))
carParameter(value::AbstractMatrix{<:Number}, name::String, unit::String) = carParameter{Array{Float64,2}}(Array{Float64,2}(value), name, unit, size(value))

# keep constructor for JuMP NonlinearExpr if you need it explicitly
# carParameter(value::NonlinearExpr, name::String, unit::String) = carParameter{NonlinearExpr}(value, name, unit, (1,1))

# Generic fallback constructor: accepts other types (VariableRef, AffExpr, custom exprs)
carParameter(value::T, name::String, unit::String) where {T} =
    begin
        sz = isa(value, AbstractVector) ? (length(value),1) :
             isa(value, AbstractMatrix) ? size(value) :
             (1,1)
        carParameter{T}(value, name, unit, sz)
    end


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
end