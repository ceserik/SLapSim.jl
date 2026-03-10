using Revise
using SLapSim
using JuMP
using Infiltrator
using Interpolations
using OrdinaryDiffEq
using DiffEqCallbacks

function timeSimulation(car::Car, result, track)
    timeVector = result.states[1:end-1, 6] #this has to be compatible with different car models will cause issues in future
    x0 = result.states[1, 1:4]
    n = result.states[1,5]
    carX = track.x[1] .- n .* sin.(track.theta[1])
    carY = track.y[1] .+ n .* cos.(track.theta[1])

    x0 = [x0; carX; carY]
    tspan = [timeVector[1], timeVector[end]]
    p = Vector{Any}(undef, 4)


    p[1] = car
    p[2] = track
    p[3] = result.controls
    p[4] = timeVector

    prob = ODEProblem(carODE_globalFrame, x0, tspan, p)

    sol = OrdinaryDiffEq.solve(prob, Tsit5(),tstops=timeVector)
    return sol
end

function timeSimulation_interpolated(car::Car, result, track)
    time(s) = result.states(s)[6] #this has to be compatible with different car models will cause issues in future
    timeVector = time.(result.path)
    time2s = linear_interpolation(timeVector,result.path) 
    initial_s = result.path[1]
    terminal_s= result.path[end]
    x0 = result.states(initial_s)[1:4]
    n_initial = result.states(initial_s)[5]

    carX = track.x[1] .- n_initial .* sin.(track.theta[1])
    carY = track.y[1] .+ n_initial .* cos.(track.theta[1])

    x0 = [x0; carX; carY]
    tspan = [time(initial_s), time(terminal_s)]
    p = Vector{Any}(undef, 4)


    p[1] = car
    p[2] = track
    p[3] = result
    p[4] = time2s

    labels = ["Vx", "Vy", "ψ", "ψ̇", "X", "Y"]
    println(rpad("t", 8), join(rpad.(labels, 10)))
    println("-"^75)
    cb = FunctionCallingCallback(; funcat = timeVector) do u, t, integrator
        vals = join(rpad.(round.(u, digits=4), 10))
        println("$(rpad(round(t,digits=3), 8))$vals")
    end

    prob = ODEProblem(carODE_globalFrame, x0, tspan, p)
    #@infiltrate
    sol = OrdinaryDiffEq.solve(prob, Tsit5(), tstops = timeVector, saveat = timeVector, callback = cb)
    return sol
end

function carODE_globalFrame(du, x, p, t)
    car = p[1]
    track = p[2]
    U = p[3]
    time2s = p[4]
    #@infiltrate
    
    u = U.controls(time2s(t))
    car.controlMapping(car, u)
    car.stateMapping(car, x)

    dx = car.carFunction(car, track, 0, nothing)

    ψ = x[3]
    rot = [cos(ψ) -sin(ψ);
        sin(ψ) cos(ψ)]

    globalFrameSpeedd = rot * [x[1]; x[2]]
    dx = [dx; globalFrameSpeedd]
    du .= dx
    return dx
end

#sol = timeSimulation(car, result, track)

function plotODESolution(sol)
    t = sol.t
    X = hcat(sol.u...)'
    labels = ["Vx", "Vy", "ψ", "ψ̇", "X", "Y"]
    println("\n  t       ", join(rpad.(labels, 10)))
    println("-"^75)
    for i in 1:length(t)
        vals = join(rpad.(round.(X[i,:], digits=4), 10))
        println("$(rpad(round(t[i],digits=3), 8))$vals")
    end
end

# Plot last two columns (X and Y coordinates) directly
#lines(getindex.(sol.u, 5), getindex.(sol.u, 6))


function plotCarPath(track::Track, result, axis = nothing)
    # take car.s and car.n and match it to track. and then add track.n
    if axis == nothing
        fig = Figure()
        axis = Axis(fig[1,1], aspect = DataAspect())
    end
    n = result.states[:, 5]
    carX = zeros(length(result.path))
    carY = zeros(length(result.path))
    for (i,s) in enumerate(result.path)
        carX[i] = track.fcurve(s)[3] .- n[i] .* sin.(track.fcurve(s)[2])
        carY[i] = track.fcurve(s)[4] .+ n[i] .* cos.(track.fcurve(s)[2])

    end

    println(axis)
    plotTrack(track, b_plotStartEnd =false, ax = axis)
    println(axis)
    lines!(axis,carX, carY, label = "Optimised")
    #return ax
end


function plotCarPath_interpolated(track::Track, result, axis = nothing)
    # take car.s and car.n and match it to track. and then add track.n
    if axis == nothing
        fig = Figure()
        axis = Axis(fig[1,1], aspect = DataAspect())
    end
    n(s) = result.states(s)[5]
    carX = zeros(length(result.path))
    carY = zeros(length(result.path))
    for (i,s) in enumerate(result.path)
        carX[i] = track.fcurve(s)[3] .- n(s) .* sin.(track.fcurve(s)[2])
        carY[i] = track.fcurve(s)[4] .+ n(s) .* cos.(track.fcurve(s)[2])

    end
#    @infiltrate
    println(axis)
    plotTrack(track, b_plotStartEnd =false, ax = axis)
    println(axis)
    lines!(axis,carX, carY, label = "Optimised")
    #return ax
end




function plotCarStates(result)
    x = result.states
    u = result.controls
    s = result.path
    fig3 = Figure()
    display(GLMakie.Screen(),fig3)
    ax = fig3[1, 1] = Axis(fig3)
    labels = ["Vx", "Vy", "ψ", "ψ̇", "n", "t"] 

    for index = 1:length(x[1,:])
        state = x[:,index]
        lines!(ax,s,state,label = labels[index], linewidth = 5)
    end
    axislegend(ax, "State Variables", position = :rt)

    fig3
    display(fig3)
    fig2 = Figure()
    display(GLMakie.Screen(),fig2)
    ax2 = [Axis(fig2[i, 1]) for i in 1:3]
    controls = ["MomentFront", "MomentRear", "Steering"] 
    for index = 1:length(u[1,:])
        control = u[:,index]
        lines!(ax2[index], s[1:end-1], control, label = controls[index], linewidth = 5)
        axislegend(ax2[index], position = :rt)
    end
    display(fig2)
end


function plotCarStates2(result)
    X = result.states
    # Plotting the results
    fig = Figure(layout = GridLayout(6, 1))
    
    # Create subplots for each variable

    labels = ["vx", "vy", "psi", "dpsi", "n", "t"]
    for (i, label) in enumerate(labels)
        ax = Axis(fig[i, 1], ylabel = label)
        lines!(ax, value.(X[:,i]), label=label)
    end
    display(GLMakie.Screen(), fig)  # This creates a new window
end