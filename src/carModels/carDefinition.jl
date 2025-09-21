
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

    # Inner constructor
   # function Car(carFunction, carParameters, controlMapping=nothing, stateMapping=nothing, mapping=nothing)
   #     new(carFunction, carParameters, controlMapping, stateMapping, mapping)
   # end
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

