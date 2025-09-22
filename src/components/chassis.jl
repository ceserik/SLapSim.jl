mutable struct Chassis
    mass::carParameter{carVar}
    hitbox::Function
end

function createCTU25chassis()
    mass = carParameter{carVar}(180.0,"mass","kg")

    function hitbox(n::carVar,track::Union{Track,Nothing},model::Union{JuMP.Model,Nothing})
        if !isnothing(model)
            #@constraint(model,n .>= -1.5)
            #@constraint(model,n .<= 1.5)
        end
    end

    chassis = Chassis(
        mass,
        hitbox
    )
    return chassis
end
