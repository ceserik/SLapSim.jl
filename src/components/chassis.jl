using Revise
mutable struct Chassis{F1}
    mass::carParameter{carVar}
    hitbox::F1
    wheelbase::carParameter{carVar}
    track::carParameter{carVar}
end

function createCTU25chassis()
    mass = carParameter{carVar}(180.0,"mass","kg")
    wheelbase = carParameter{carVar}(1.525,"wheelbase","m")
    track = carParameter{carVar}(1.2,"track","m")

    function hitbox(n::carVar,track::Union{Track,Nothing},model=nothing)
        greaterContraint(n, -track.widthR[1], model)
        lessContraint(n, track.widthL[1], model)
    end

    chassis = Chassis(
        mass,
        hitbox,
        wheelbase,
        track
    )
    return chassis
end
