
using Interpolations

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

