
function lessContraint(a::carVar, b::carVar, model::Union{JuMP.Model,Nothing}=nothing)
    if isnothing(model)
        return min(a, b)
    else
        @constraint(model, a <= b)
        return a
    end
end

function greaterContraint(a::carVar, b::carVar, model::Union{JuMP.Model,Nothing}=nothing)
    if isnothing(model)
        return max(a, b)
    else
        @constraint(model, a >= b)
        return a
    end
end

function equalConstraint(a::carVar, b::carVar, model::Union{JuMP.Model,Nothing}=nothing)
    if isnothing(model)
        return b
    else
        @constraint(model, a == b)
        return a
    end
end




mutable struct Car2
    drivetrain
    aero
    suspension
    chassis
    wheelAssemblies
    carParameters
    carODE
    controlMapping::Function
    stateMapping::Function
end

mutable struct Drivetrain
    motors::Vector{Motor}
    gearboxes::Vector{Gearbox}
    tires::Vector{Tire}
    accumulators::Accumulator
end

# Define the Car struct with parameters
mutable struct Car
    carFunction::Function
    carParameters::CarParameters
    controlMapping::Function
    stateMapping::Function
    mapping::Function
    drivetrain::Drivetrain
    aero::Aero
    suspension::Suspension
    chassis::Chassis
    wheelAssemblies::Vector{WheelAssembly}
    state_desc::Union{Vector{VarEntry}, Nothing}
    control_desc::Union{Vector{VarEntry}, Nothing}
end






Base.show(io::IO, ::MIME"text/plain", obj::Car) = prettyPrintComponent(io, obj)
Base.show(io::IO, ::MIME"text/plain", obj::Drivetrain) = prettyPrintComponent(io, obj)

"""
Generic pretty printing for car components.
Usage: import this function and add:
    Base.show(io::IO, ::MIME"text/plain", obj::YourType) = prettyPrintComponent(io, obj)
"""
function prettyPrintComponent(io::IO, obj; indent::Int=0)
    # List of types that should get pretty printing
    pretty_print_types = [Tire, Motor, Gearbox, Chassis, Drivetrain, Car, Car2, WheelAssembly, CarParameters]

    # Handle types that shouldn't get pretty printing
    if !any(T -> obj isa T, pretty_print_types)
        print(io, obj)
        return
    end

    # Handle car components with pretty printing
    T = typeof(obj)
    prefix = "  " ^ indent
    inner = "  " ^ (indent + 1)
    println(io, "$(T)(")

    if fieldcount(T) > 0
        max_width = maximum(length(string(field)) for field in fieldnames(T))
        for field in fieldnames(T)
            value = getfield(obj, field)
            value isa Function && continue
            field_str = rpad(string(field), max_width)
            if any(PT -> value isa PT, pretty_print_types)
                print(io, "$(inner)$(field_str) = ")
                prettyPrintComponent(io, value; indent=indent+1)
                println(io)
            else
                println(io, "$(inner)$(field_str) = ", value)
            end
        end
    end

    print(io, "$(prefix))")
end

