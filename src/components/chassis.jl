using Revise
mutable struct Chassis{F1}
    mass::carParameter{carVar}
    hitbox::F1
    wheelbase::carParameter{carVar}
    track::carParameter{carVar}
    CoG_X_pos::carParameter{carVar}
    CoG_Y_pos::carParameter{carVar}
end

function createCTU25chassis()
    mass = carParameter{carVar}(280.0,"mass","kg")
    wheelbase = carParameter{carVar}(1.525,"wheelbase","m")
    track = carParameter{carVar}(1.2,"track","m")
    CoG_X_pos = carParameter{carVar}(0.5,"ratio of CoG on front","m")
    CoG_Y_pos = carParameter{carVar}(0.5,"ratio of CoG on left","m")
    function hitbox(n::carVar,track::Union{Track,Nothing},model=nothing)
        greaterContraint(n, -track.widthR[1]+0.6, model)
        lessContraint(n, track.widthL[1]-0.6, model)
    end

    chassis = Chassis(
        mass,
        hitbox,
        wheelbase,
        track,
        CoG_X_pos,
        CoG_Y_pos
    )
    return chassis
end

"""
    draw!(ax, chassis::Chassis, x, y, ψ)

Draw the chassis body as a rectangle centered at (x,y) with heading ψ.
"""
function draw!(ax, chassis::Chassis, x, y, ψ)
    wb = chassis.wheelbase.value
    tw = chassis.track.value
    c = _rect_corners(x, y, wb + 0.2, tw + 0.1, ψ)
    poly!(ax, Point2f.(eachrow(c)); color=(:steelblue, 0.3), strokecolor=:steelblue, strokewidth=2)
end
