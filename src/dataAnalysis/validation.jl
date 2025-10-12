using Revise
using SLapSim
using JuMP
using Infiltrator
using Interpolations
using OrdinaryDiffEq

function timeSimulation(car::Car,result::SLapSim.Result,track)
    timeVector = result.states[1:end-1,6] #this has to be compatible with different car models
    x0 = result.states[1,1:4]
    x0 =[x0;-11;-20]
    tspan = [timeVector[1],timeVector[end]]
    p = Vector{Any}(undef, 4)
    
    
    p[1] = car
    p[2] = track
    p[3] = result.controls
    p[4] = timeVector

    prob = ODEProblem(carODE_globalFrame,x0,tspan,p)
    sol = solve(prob, Tsit5())
    return sol
end

function carODE_globalFrame(du,x,p,t)
    car = p[1]
    track = p[2]
    U = p[3]
    timeVector = p[4]
    interp_linear = linear_interpolation((timeVector,1:3),U)
    u = interp_linear(t,1:3)
    @infiltrate
    car.controlMapping(car,u)
    car.stateMapping(car,x)
    
    dx = car.carFunction(car,track,0,nothing)

    ψ = x[3]
    rot = [ cos(ψ) -sin(ψ);
            sin(ψ) cos(ψ)]

    globalFrameSpeedd = rot* [x[1];x[2]]
    dx = [dx;globalFrameSpeedd]
    du .=dx
    return dx
end

sol = timeSimulation(car,result,track)

# Plot last two columns (X and Y coordinates) directly
lines(getindex.(sol.u, 5), getindex.(sol.u, 6))