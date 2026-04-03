mutable struct Aero{F1,F2,F3}
    CL::carParameter{carVar}
    CD::carParameter{carVar}
    CoP::carParameter{carVar}
    compute::F1
    setupObservables::F2
    updateObservables::F3
end
#Base.show(io::IO, ::MIME"text/plain", obj) = prettyPrintComponent(io, obj)

const RHO_SEA_LEVEL = 1.225

function createBasicAero(CL_a::Float64 = -5.0, CD_a::Float64 = 2.0, CoP_a::Float64 = 0.5)
    CL = carParameter{carVar}(CL_a,"Lift coeffcient","-")
    CD = carParameter{carVar}(CD_a,"Drag coeffcient","-")
    CoP = carParameter{carVar}(CoP_a,"Cener of pressure on front","-")
    function compute(vx::carVar, rho::Float64=RHO_SEA_LEVEL)
        downforce = -0.5 * rho * CL.value * vx^2
        drag = -0.5 * rho * CD.value * vx^2
        return (downforce=downforce, drag=drag)
    end
    function setupObservables(ax::Axis)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        front_obs = Observable(dummy)
        rear_obs  = Observable(dummy)
        poly!(ax, front_obs; color=:transparent, strokecolor=:transparent, strokewidth=1.5)
        poly!(ax, rear_obs;  color=:transparent, strokecolor=:transparent, strokewidth=1.5)
        return (front=front_obs, rear=rear_obs)
    end

    function updateObservables(obs::NamedTuple, x::Float64, y::Float64, ψ::Float64, wheelbase::Float64, trackwidth::Float64)
        R = _rotmat2d(ψ)
        fc = R * [wheelbase / 2 + 0.15; 0.0] .+ [x, y]
        obs.front[] = _rect_points(fc[1], fc[2], 0.05, trackwidth + 0.3, ψ)
        rc = R * [-wheelbase / 2 - 0.15; 0.0] .+ [x, y]
        obs.rear[] = _rect_points(rc[1], rc[2], 0.05, trackwidth + 0.2, ψ)
    end

    aero = Aero(
        CL,
        CD,
        CoP,
        compute,
        setupObservables,
        updateObservables,
    )
    return aero
end
