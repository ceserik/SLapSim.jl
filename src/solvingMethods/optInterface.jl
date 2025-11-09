
using Revise
using Infiltrator
using SLapSim
using JuMP
using OrdinaryDiffEq
using DifferentiationInterface

struct Result
    states::Matrix{Float64}
    controls::Matrix{Float64}
    path::Vector{Float64}
end


function initializeSolution(car::Car,track::Track,sampleDistances::Vector{Float64})
    span2 = [sampleDistances[1],sampleDistances[end]]
    x0 = [2.0, 0.0, track.theta[1], 0.0, 0, 0.0]
    steeringP = 10
    velocityP = 3
    vref = 2

    prob = ODEProblem(carODE_path2, x0, span2,[car,track,steeringP,velocityP,vref]);
    sol = solve(prob,Tsit5(),saveat=sampleDistances,reltol=1e-5, abstol=1e-5);

    x = hcat(sol.u...)'
    s = sol.t
    steering = x[:,5] .* steeringP
    torque = (vref .- x[:,1]) * velocityP

    u = zeros(length(sampleDistances),Int64(car.carParameters.nControls.value))
    u[:,1] = torque
    u[:,3] = steering

    initialization = Result(
        x,
        u[1:end-1,:],
        s
    )
    return initialization
end;
fig_interp = Figure()


function carODE_path(
    car::Car,track::Track,
    k::Union{Int64,Float64},
    u::Union{Vector{VariableRef},
    Vector{Float64}},
    x::Union{Vector{VariableRef},Vector{Float64}},
    model::Union{JuMP.Model,Nothing},
    make_constraints::Bool)
    #car.mapping(car,u,x)
    car.controlMapping(car,u)
    car.stateMapping(car,x)
    dzds = time2path(car,track,k,model,make_constraints) #time2path(s,instantCarParams,track,car)
    return dzds
end;

#poradie car,track,k,s,u,x


function carODE_path2(du, x, p, s)
    car = p[1];
    track = p[2]
    steeringP = p[3]
    velocityP = p[4]
    vref = p[5]

    control = [(vref - car.carParameters.velocity.value[1])*velocityP, 0.0, x[5]*steeringP]

    dx = carODE_path(car,track,s, control, x,nothing,false)  ; # carF must return a vector matching length(x)
    du .= dx;
end;


#initializeSolution(1,2)
## here I need to define transfomation of ODE with respect to time to ODE with respect to path
function time2path(car::Car,track::Track,k::Union{Int64,Float64},model::Union{Nothing,JuMP.Model},make_constraints::Bool=false)
    #track.mapping(track,instantTrack,s)
    dxdt = car.carFunction(car,track,k,model,make_constraints)
    v_x = car.carParameters.velocity.value[1]
    v_y = car.carParameters.velocity.value[2]
    psi = car.carParameters.psi.value
    n   = car.carParameters.n.value

    th = track.fcurve(k)[2]
    C = track.fcurve(k)[1]
    
    epsilon = psi - th

    Sf = (1 - n*C) ./ (v_x.*cos(epsilon) - v_y.*sin(epsilon));
    dndt = v_x.*sin(epsilon) + v_y.*cos(epsilon);

    dzds = [
    Sf .* dxdt[1];    #dvx/ds
    Sf .* dxdt[2];    #dvy/ds
    Sf .* dxdt[3];    #dpsi/ds
    Sf .* dxdt[4];    #ddpsi/ds
    Sf .* dndt;       #dn/ds
    Sf;               #dt/ds 
    # this should be capable of accepting more states
    ]

    #Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
    return dzds
end




function findOptimalTrajectory(track::Track,car::Car,model::JuMP.Model,sampleDistances, initialization=nothing)
    N = length(sampleDistances)

    #fill track parameters which are constant along track, should be automatized
    track.rho    = fill(track.rho[1],N)
    track.μ      = fill(track.μ[1],N)

    nStates =  Int(car.carParameters.nStates.value)
    s = sampleDistances
    #create inputs for car model #create instance of track parameters

    #determine sizes of inputs and states
    nControls = Int(round(car.carParameters.nControls.value))
    
    function f(x,u,s,make_constraints)
        dxds = carODE_path(car,track,s,u,x,model,make_constraints)
        return dxds
    end

    stages = 2
    lobotom= createLobattoIIIA(stages,f)
    xd = lobotom.createConstraints(f,6,3,track.fcurve,s,model,initialization.states,initialization.controls)
    X     = xd[2]
    U     = xd[3][1:end-1,:]
    s_all = xd[6]
    @objective(model,Min,X[end,6])


    @constraint(model,X[1:end,1] .>= 0) #vx
    #@constraint(model,X[1,1] .== 5) # intial vx
    @constraint(model,X[1,2] .== 0) # intial vy
    @constraint(model,X[1,3] .== track.theta[1]) # intial heading
    #@constraint(model,X[end,3].== track.theta[end])
    #@constraint(model,X[end,5].== 0)
    @constraint(model,X[1,6] .>= 0) # final time
    @constraint(model,diff(X[:,6]) .>=0) #time goes forward


    
    optimize!(model)

    out = Result(value.(X),value.(U),s_all)
    return out
end

function simulateInTime(car::Car,track::Track,result::Result)





end