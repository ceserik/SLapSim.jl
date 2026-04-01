mutable struct Aero{F1}
    CL::carParameter{carVar}
    CD::carParameter{carVar}
    CoP::carParameter{carVar}
    compute::F1
end
#Base.show(io::IO, ::MIME"text/plain", obj) = prettyPrintComponent(io, obj)

const RHO_SEA_LEVEL = 1.225

function createBasicAero()
    CL = carParameter{carVar}(-5.0,"Lift coeffcient","-")
    CD = carParameter{carVar}(2.0,"Drag coeffcient","-")
    CoP = carParameter{carVar}(0.5,"Cener of pressure on front","-")
    function compute(vx, rho=RHO_SEA_LEVEL)
        downforce = 0.5 * rho * CL.value * vx^2
        drag = -0.5 * rho * CD.value * vx^2
        return (downforce=downforce, drag=drag)
    end
    aero = Aero(
        CL,
        CD,
        CoP,
        compute
    )
    return aero

end


function setup_observables!(ax, aero::Aero)
    dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
    front_obs = Observable(dummy)
    rear_obs  = Observable(dummy)
    poly!(ax, front_obs; color=:transparent, strokecolor=:transparent, strokewidth=1.5)
    poly!(ax, rear_obs;  color=:transparent, strokecolor=:transparent, strokewidth=1.5)
    return (front=front_obs, rear=rear_obs)
end

function update_observables!(obs, aero::Aero, x, y, ψ, wheelbase, trackwidth)
    R = _rotmat2d(ψ)
    fc = R * [wheelbase / 2 + 0.15; 0.0] .+ [x, y]
    obs.front[] = _rect_points(fc[1], fc[2], 0.05, trackwidth + 0.3, ψ)
    rc = R * [-wheelbase / 2 - 0.15; 0.0] .+ [x, y]
    obs.rear[] = _rect_points(rc[1], rc[2], 0.05, trackwidth + 0.2, ψ)
end
