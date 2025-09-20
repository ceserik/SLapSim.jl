
function lessContraint(a,b,model=nothing)
    if isnothing(model)
        a = min(a,b)
    else
        @constraint(model,a<=b)
    end
    return a
end

function greaterContraint(a,b,model=nothing)
    if isnothing(model)
        a = max(a,b)
    else
        @constraint(model,a>=b)
    end
    return a
end

function equalConstraint(a,b)
    if isnothing(model)
        a = b
    else
        @constraint(model,a==b)
    end    
    return a
end


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


mutable struct carParameters
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

