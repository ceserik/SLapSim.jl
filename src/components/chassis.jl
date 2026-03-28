using Revise
mutable struct Chassis
    mass::carParameter{carVar}
    hitbox::Function
    wheelbase::carParameter{carVar}
    track::carParameter{carVar}
end

function createCTU25chassis()
    mass = carParameter{carVar}(180.0,"mass","kg")
    wheelbase = carParameter{carVar}(1.525,"wheelbase","m")
    track = carParameter{carVar}(1.2,"track","m")

    function hitbox(n::carVar,track::Union{Track,Nothing},model::Union{JuMP.Model,Nothing})
        if !isnothing(model)
            @constraint(model,n >= -track.widthR[1])
            @constraint(model,n <= track.widthL[1])
        end
    end

    chassis = Chassis(
        mass,
        hitbox,
        wheelbase,
        track
    )
    return chassis
end
