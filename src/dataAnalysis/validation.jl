using SLapSim
using JuMP

using Interpolations
using OrdinaryDiffEq
using DiffEqCallbacks
using LinearAlgebra
using LaTeXStrings



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
    header = rpad("t", 8) * join(rpad.(labels, 10))
    println(header)
    println("-"^length(header))
    cb = FunctionCallingCallback(; funcat = t_unique) do u, t, integrator
        vals = join(rpad.(round.(u, digits=4), 10))
        print("\r$(rpad(round(t,digits=3), 8))$vals")
    end

    prob = ODEProblem(carODE_globalFrame, x0, tspan, p)
    sol = OrdinaryDiffEq.solve(prob, AutoTsit5(Rodas4(autodiff = AutoFiniteDiff())),  saveat=t_unique, callback=cb, reltol=1e-3, abstol=1e-3)

    return sol
end

function carODE_globalFrame(du, x, p, t)
    car = p[1]
    track = p[2]
    U = p[3]
    time2s = p[4]
    #@infiltrate
    
    u = U.controls(time2s(t))
    car.controlMapping(u)
    car.stateMapping(x)

    dx = car.carFunction(track, nothing)

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
    s = result.path[1]:stepsize:result.path[end]
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

"""
    getSegmentErrors(problem) -> Vector{Float64}

One scaled error per mesh segment. Integrates the ODE across each segment using
the optimized polynomial controls, then compares the simulated state at the
right endpoint to the optimizer's state there. States are normalized by
`get_scales(state_descriptor)` so each contributes comparably.
"""
function getSegmentErrors(problem; method::Symbol = :ode)
    edges = problem.optiResult.segment_edges
    nseg  = length(edges) - 1
    errs  = zeros(Float64, nseg)
    x_scale = get_scales(problem.car.carParameters.state_descriptor)
    for i = 1:nseg
        print("\rSegment error ($method): $i/$nseg")
        seg = [edges[i], edges[i+1]]
        errs[i] = method === :defect ?
            getError_defect(seg, problem) :
            getError(seg, problem; x_scale=x_scale)
    end
    println()
    return errs
end


"""ODE-integration error between every pair of stored nodes; plot vs s."""
function getNodeErrors(problem; iteration::Int=0, savepath::Union{Nothing,String}=nothing)
    p = problem.optiResult.path
    xs = get_scales(problem.car.carParameters.state_descriptor)
    errs = [p[i] == p[i+1] ? 0.0 : getError([p[i], p[i+1]], problem; x_scale=xs) for i in 1:length(p)-1]
    s_mid = [(p[i] + p[i+1]) / 2 for i in 1:length(p)-1]
    fig = scatter(s_mid, errs; axis=(; yscale=log10, title="Node-interval errors — iter $iteration", xlabel="s", ylabel="Node Error"))
    if savepath !== nothing
        save(savepath, fig)
    end
    display(GLMakie.Screen(), fig)
    return errs, s_mid
end



function getError(s, problem; x_scale=nothing)
    x0 = problem.optiResult.states(s[1])
    ode(du, x, p, s) = du .= carODE_path(p[1], p[2], s, p[3](s), x, nothing)
    prob = ODEProblem(ode, x0, (s[1], s[2]), (problem.car, problem.track, problem.optiResult.controls))

    # Loose tolerance + step cap: refinement only needs a coarse error estimate, and a
    # pathological polynomial (high order, long segment) can otherwise stall Rodas4.
    sol = OrdinaryDiffEq.solve(prob, Rodas4(autodiff=AutoFiniteDiff());
                               reltol=1e-3, abstol=1e-3, maxiters=1_000, verbose=false)
    if sol.retcode != ReturnCode.Success
        return Inf  # treat as "needs refinement"
    end
    final_states_time_sim  = sol.u[end]
    final_states_optimized = problem.optiResult.states(s[2])
    diff = final_states_time_sim .- final_states_optimized
    if x_scale !== nothing
        diff = diff ./ x_scale
    end
    return norm(diff)
end

# Vector-valued Romberg quadrature on [a,b]: trapezoidal rule with successive
# halving + Richardson extrapolation. Returns ∫_a^b f(s) ds componentwise.
# Keeps only the previous and current row of the Romberg tableau.
function romberg_vec(f, a::Float64, b::Float64; tol::Float64 = 1e-3, maxiter::Int = 6)
    h    = b - a
    fa   = f(a)
    fb   = f(b)
    dim  = length(fa)
    Rprev = [0.5 * h .* (fa .+ fb)]
    Rcur  = Vector{Vector{Float64}}(undef, 0)
    for k in 2:maxiter
        h /= 2
        n = 1 << (k - 2)              # 2^(k-2)
        s_acc = zeros(dim)
        @inbounds for ii in 1:n
            s_acc .+= f(a + (2ii - 1) * h)
        end
        Rcur = Vector{Vector{Float64}}(undef, k)
        Rcur[1] = 0.5 .* Rprev[1] .+ h .* s_acc
        @inbounds for j in 2:k
            pow = 1 << (2 * (j - 1))  # 4^(j-1)
            Rcur[j] = (pow .* Rcur[j-1] .- Rprev[j-1]) ./ (pow - 1)
        end
        if maximum(abs.(Rcur[k] .- Rprev[k-1])) < tol * maximum(abs.(Rcur[k]) .+ 1.0)
            return Rcur[k]
        end
        Rprev = Rcur
    end
    return Rcur[end]
end

# Per Betts §4.7, eq. (4.154)–(4.156): integrate |ε_i(s)| = |ỹ̇_i − f_i| over
# the segment via Romberg quadrature, return the max relative error
# max_i η_i / (w_i + 1) using per-state scale weights.
#
# Defaults tuned for mesh-refinement classification (not machine precision).
# Tighten tol/maxiter if you want a true Betts-grade integral.
function getError_defect(s_segment, problem; tol::Float64 = 1e-3, maxiter::Int = 6)
    a, b = Float64(s_segment[1]), Float64(s_segment[2])
    x_scale     = get_scales(problem.car.carParameters.state_descriptor)
    states      = problem.optiResult.states
    controls    = problem.optiResult.controls
    car         = problem.car
    track       = problem.track
    state_deriv = problem.optiResult.state_deriv

    h_fd = 1e-5 * (b - a)
    inv_2h = 1 / (2h_fd)

    if state_deriv === nothing
        eps_abs = function (s)
            x    = states(s)
            u    = controls(s)
            f    = carODE_path(car, track, s, u, x, nothing)
            dxds = (states(s + h_fd) .- states(s - h_fd)) .* inv_2h
            return abs.(dxds .- f)
        end
    else
        eps_abs = function (s)
            x    = states(s)
            u    = controls(s)
            f    = carODE_path(car, track, s, u, x, nothing)
            return abs.(state_deriv(s) .- f)
        end
    end

    eta = romberg_vec(eps_abs, a, b; tol = tol, maxiter = maxiter)
    return maximum(eta ./ (x_scale .+ 1))
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

function plotErrorsOnTrack2D(problem; axis=nothing, colormap=:turbo, itp=nothing)
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


function plot_on_path(problem,itp,legend; axis=nothing, colormap=:turbo)
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
    plt = lines!(axis, carX, carY; color = vis_vars, colormap = colormap, linewidth = 10)
    if created
        Colorbar(fig[1,2], plt; label = legend, width = 25)
        display(GLMakie.Screen(), fig)
        return fig, axis, plt
    else
        return axis, plt
    end
end


function plot_controls_on_path(problem, optiResult_interp; error_itp=nothing)
    n_controls = Int(problem.car.carParameters.nControls.value)
    desc = problem.car.carParameters.control_descriptor
    plots = [(s -> optiResult_interp.controls(s)[i],
              desc !== nothing ? desc[i].name : "control_$i") for i in 1:n_controls]
    if !isnothing(error_itp)
        push!(plots, (s -> error_itp(s), "error"))
    end

    n = length(plots)
    fig = Figure(size=(900, 450 * ceil(Int, n / 2)))
    for i in 1:n
        row   = div(i - 1, 2) + 1
        col   = mod(i - 1, 2) + 1
        ax    = Axis(fig[row, 2*col-1], aspect=DataAspect())
        _, plt = plot_on_path(problem, plots[i][1], plots[i][2]; axis=ax)
        Colorbar(fig[row, 2*col], plt; label=plots[i][2], width=25)
    end
    display(GLMakie.Screen(), fig)
    return fig
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