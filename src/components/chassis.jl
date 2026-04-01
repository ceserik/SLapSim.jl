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
    function hitbox(n::carVar,racetrack::Union{Track,Nothing},model=nothing)
        greaterContraint(n, -racetrack.widthR[1] + track.value / 2, model)
        lessContraint(n, racetrack.widthL[1] - track.value / 2, model)
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



function setup_observables!(ax, chassis::Chassis)
    dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
    obs = Observable(dummy)
    poly!(ax, obs; color=:transparent, strokecolor=:transparent, strokewidth=0.5)
    return obs
end

function update_observables!(obs, chassis::Chassis, x, y, ψ)
    wb = chassis.wheelbase.value
    tw = chassis.track.value
    obs[] = _rect_points(x, y, wb + 0.2, tw + 0.1, ψ)
end
