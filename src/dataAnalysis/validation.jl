using Revise
using SLapSim
using JuMP
using Infiltrator
using Interpolations
using OrdinaryDiffEq
using DiffEqCallbacks
using LinearAlgebra
using LaTeXStrings

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
    timeVector = accumulate(max, time.(result.path))  # enforce monotonicity against optimizer tolerance violations
    # remove duplicates (from optimizer tolerance) so interpolation knots are strictly increasing
    unique_mask = [true; diff(timeVector) .> 0]
    t_unique   = timeVector[unique_mask]
    s_unique   = result.path[unique_mask]
    # Flat() extrapolation prevents BoundsError on RK sub-steps past tspan[end]
    #i should find cleaner way of fixing this
    time2s = extrapolate(interpolate((t_unique,), s_unique, Gridded(Linear())), Flat())

    initial_s = result.path[1]
    terminal_s = result.path[end]
    x0 = result.states(initial_s)[1:4]
    n_initial = result.states(initial_s)[5]

    carX = track.x[1] .- n_initial .* sin.(track.theta[1])
    carY = track.y[1] .+ n_initial .* cos.(track.theta[1])

    x0 = [x0; carX; carY]
    tspan = [t_unique[1], t_unique[end]]
    p = Vector{Any}(undef, 4)

    p[1] = car
    p[2] = track
    p[3] = result
    p[4] = time2s

    labels = ["Vx", "Vy", "ψ", "ψ̇", "X", "Y"]
    println(rpad("t", 8), join(rpad.(labels, 10)))
    println("-"^75)
    cb = FunctionCallingCallback(; funcat = t_unique) do u, t, integrator
        vals = join(rpad.(round.(u, digits=4), 10))
        println("$(rpad(round(t,digits=3), 8))$vals")
    end

    prob = ODEProblem(carODE_globalFrame, x0, tspan, p)
    sol = OrdinaryDiffEq.solve(prob, Tsit5(), tstops = t_unique, saveat = t_unique, callback = cb)
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
    created = false
    if axis === nothing
        fig = Figure()
        axis = Axis(fig[1,1], aspect = DataAspect())
        created = true
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
    
    if created
        display(GLMakie.Screen(), fig)
    end
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




function plotCarStates_interp(result, stepsize)
    s = collect(result.path[1]:stepsize:result.path[end])
    #s  = result.path
    x = hcat(result.states.(s)...)'   # n_samples × nStates
    u = hcat(result.controls.(s)...)'  # n_samples × nControls

    fig3 = Figure()
    display(GLMakie.Screen(), fig3)
    labels = ["Vx", "Vy", "ψ", "ψ̇", "n", "t"]
    ax3 = [Axis(fig3[i, 1], ylabel = labels[i]) for i in 1:size(x, 2)]

    for index = 1:size(x, 2)
        lines!(ax3[index], s, x[:, index], label = labels[index], linewidth = 5)
        axislegend(ax3[index], position = :rt)
    end
    display(fig3)

    fig2 = Figure()
    display(GLMakie.Screen(), fig2)
    ax2 = [Axis(fig2[i, 1]) for i in 1:size(u, 2)]
    ctrl_labels = ["MomentFront", "MomentRear", "Steering"]
    for index = 1:size(u, 2)
        lines!(ax2[index], s, u[:, index], label = ctrl_labels[index], linewidth = 5)
        axislegend(ax2[index], position = :rt)
    end
    display(fig2)
end


#This function applies states and control from optimisation result and using Tsit5 predicts next states, the difference 
#between states from optimisation and states from Tsit5 simulation is compared and used as output.

function getErrors(problem)
    s = problem.optiResult.path
    n = length(s)
    error_vector = zeros(Float64, n)
    for i = 1:n-1
        
        error = getError([s[i], s[i+1]],problem)
        error_vector[i] = error
        #println(error)
    end
    error_vector[end] = error_vector[end-1]
    itp = linear_interpolation(s, error_vector)
    return itp
end



function getError(s, problem)
    x0 = problem.optiResult.states(s[1])
    ode(du, x, p, s) = du .= carODE_path(p[1], p[2], s, p[3](s), x, nothing)
    prob = ODEProblem(ode, x0, (s[1], s[2]), (problem.car, problem.track, problem.optiResult.controls))
    
    sol = OrdinaryDiffEq.solve(prob, Rodas4(autodiff = AutoFiniteDiff()), reltol=1e-5, abstol=1e-5)
    final_states_time_sim = sol.u[end]
    final_states_optimized= problem.optiResult.states(s[2])
    #@infiltrate
    # return scalar L2 norm (single number)
    return norm(final_states_time_sim .- final_states_optimized)
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

function plotErrorsOnTrack2D(problem; axis=nothing, colormap=:viridis, itp=nothing)
    track = problem.track
    if itp === nothing
        itp = getErrors(problem)          # interpolation from getErrors
    end
    s = problem.optiResult.path
    errs = itp.(s)

    # compute car trajectory (offset from centerline) like plotCarPath_interpolated
    n(s) = problem.optiResult.states(s)[5]
    carX = zeros(length(s))
    carY = zeros(length(s))
    for (i, si) in enumerate(s)
        fc = track.fcurve(si)
        carX[i] = fc[3] .- n(si) * sin(fc[2])
        carY[i] = fc[4] .+ n(si) * cos(fc[2])
    end

    created = false
    if axis === nothing
        fig = Figure()
        axis = Axis(fig[1,1], aspect = DataAspect())
        created = true
    end
    

    plotTrack(track, b_plotStartEnd=false, ax = axis)   # draw base track
    plt = lines!(axis, carX, carY; color = errs, colormap = colormap, linewidth = 4)
    if created
        Colorbar(fig[1,2], plt; label = "error", width = 25)
        display(GLMakie.Screen(), fig)
        return fig, axis, plt
    else
        return axis, plt
    end
end


function plot_on_path(problem,itp,legend; axis=nothing, colormap=:viridis)
    track = problem.track

    #itp = getErrors(problem)              # interpolation from getErrors
    s = problem.optiResult.path
    vis_vars = itp.(s)

    # compute car trajectory (offset from centerline) like plotCarPath_interpolated
    n(s) = problem.optiResult.states(s)[5]
    carX = zeros(length(s))
    carY = zeros(length(s))
    for (i, si) in enumerate(s)
        fc = track.fcurve(si)
        carX[i] = fc[3] .- n(si) * sin(fc[2])
        carY[i] = fc[4] .+ n(si) * cos(fc[2])
    end

    created = false
    if axis === nothing
        fig = Figure()
        axis = Axis(fig[1,1], aspect = DataAspect())
        created = true
    end

    plotTrack(track, b_plotStartEnd=false, ax = axis)   # draw base track
    plt = lines!(axis, carX, carY; color = vis_vars, colormap = colormap, linewidth = 4)
    if created
        Colorbar(fig[1,2], plt; label = legend, width = 25)
        display(GLMakie.Screen(), fig)
        return fig, axis, plt
    else
        return axis, plt
    end
end


function get_sampling_density(segment_edges)
    h = diff(segment_edges)
    density = 1 ./ h
    left_edges = segment_edges[1:end-1]

    #fig = plot(left_edges, density,
    #    xlabel="Track position (m)",
    #    ylabel="Segments / m",
    #    title="Sampling density",
    #    seriestype=:steppre,
    #    label="density")
    #display(fig)

    density_itp = extrapolate(interpolate((left_edges,), density, Gridded(Linear())), Flat())

    return density_itp
end