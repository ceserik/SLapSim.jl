mutable struct Aero
    CL::carParameter{carVar}
    CD::carParameter{carVar}
    CoP::carParameter{carVar}
end
#Base.show(io::IO, ::MIME"text/plain", obj) = prettyPrintComponent(io, obj)

function createBasicAero()
    CL = carParameter{carVar}(-5.0,"Lift coeffcient","-")
    CD = carParameter{carVar}(2.0,"Drag coeffcient","-")
    CoP = carParameter{carVar}(0.5,"Cener of pressure on front","-")
    aero = Aero(
        CL,
        CD,
        CoP
    )
    return aero

end

"""
    draw!(ax, aero::Aero, x, y, ψ, wheelbase, trackwidth)

Draw front and rear wings. Needs wheelbase and trackwidth from chassis for positioning.
"""
function draw!(ax, aero::Aero, x, y, ψ, wheelbase, trackwidth)
    R = _rotmat2d(ψ)
    # Front wing
    fc = R * [wheelbase / 2 + 0.15; 0.0] .+ [x, y]
    c = _rect_corners(fc[1], fc[2], 0.05, trackwidth + 0.3, ψ)
    poly!(ax, Point2f.(eachrow(c)); color=:red, strokecolor=:darkred, strokewidth=1.5)
    # Rear wing
    rc = R * [-wheelbase / 2 - 0.15; 0.0] .+ [x, y]
    c = _rect_corners(rc[1], rc[2], 0.05, trackwidth + 0.2, ψ)
    poly!(ax, Point2f.(eachrow(c)); color=:red, strokecolor=:darkred, strokewidth=1.5)
end

function setup_observables!(ax, aero::Aero)
    dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
    front_obs = Observable(dummy)
    rear_obs  = Observable(dummy)
    poly!(ax, front_obs; color=:red, strokecolor=:darkred, strokewidth=1.5)
    poly!(ax, rear_obs;  color=:red, strokecolor=:darkred, strokewidth=1.5)
    return (front=front_obs, rear=rear_obs)
end

function update_observables!(obs, aero::Aero, x, y, ψ, wheelbase, trackwidth)
    R = _rotmat2d(ψ)
    fc = R * [wheelbase / 2 + 0.15; 0.0] .+ [x, y]
    obs.front[] = _rect_points(fc[1], fc[2], 0.05, trackwidth + 0.3, ψ)
    rc = R * [-wheelbase / 2 - 0.15; 0.0] .+ [x, y]
    obs.rear[] = _rect_points(rc[1], rc[2], 0.05, trackwidth + 0.2, ψ)
end
