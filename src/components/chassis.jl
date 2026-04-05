mutable struct Chassis{F1,F2,F3}
    mass::carParameter{carVar}
    hitbox::F1
    wheelbase::carParameter{carVar}
    track::carParameter{carVar}
    CoG_X_pos::carParameter{carVar}
    CoG_Y_pos::carParameter{carVar}
    CoG_Z_pos::carParameter{carVar}
    width::carParameter{carVar}
    setupObservables::F2
    updateObservables::F3
end

function createCTU25chassis(; mass_p::Float64 = 280.0, wheelbase_p::Float64 = 1.525, track_p::Float64 = 1.2, CoGx::Float64 = 0.5, CogY::Float64 = 0.5, width_p::Float64 = 1.525)
    mass = carParameter{carVar}(mass_p,"mass","kg")
    wheelbase = carParameter{carVar}(wheelbase_p,"wheelbase","m")
    track = carParameter{carVar}(track_p,"track","m")
    width = carParameter{carVar}(width_p,"track","m")
    CoG_X_pos = carParameter{carVar}(CoGx,"ratio of CoG on front","-",:tunable)
    CoG_Y_pos = carParameter{carVar}(CogY,"ratio of CoG on left","-",:tunable)
    CoG_Z_pos = carParameter{carVar}(0.2,"CoG height in z from road","m",:tunable)
    function hitbox(n::carVar,racetrack::Union{Track,Nothing},model=nothing)
        #greaterContraint(n, -racetrack.widthR[1] + track.value / 2, model)
        #lessContraint(n, racetrack.widthL[1] - track.value / 2, model)
    end

    function setupObservables(ax)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        obs = Observable(dummy)
        poly!(ax, obs; color=:transparent, strokecolor=:transparent, strokewidth=0.5)
        return obs
    end

    function updateObservables(obs::Observable, x::Float64, y::Float64, ψ::Float64)
        obs[] = _rect_points(x, y, wheelbase.value + 0.2, track.value + 0.1, ψ)
    end

    chassis = Chassis(
        mass,
        hitbox,
        wheelbase,
        track,
        CoG_X_pos,
        CoG_Y_pos,
        CoG_Z_pos,
        width,
        setupObservables,
        updateObservables,
    )
    return chassis
end
