
using SLapSim
using JuMP
using OrdinaryDiffEq
using DifferentiationInterface
using Interpolations


struct Result
    states::Matrix{Float64}
    controls::Matrix{Float64}
    path::Vector{Float64}
end


mutable struct Problem_config
    car::Union{Nothing, Car}
    track::Union{Nothing, Track}
    model::Union{Nothing, JuMP.Model}
    optiResult
    params
end

struct Result_interpolation
    states
    controls
    path
end

function make_result_interpolation(x::AbstractMatrix{Float64}, u::AbstractMatrix{Float64}, s_states::Vector{Float64}, s_controls::Vector{Float64})
    state_interps = [linear_interpolation(s_states, x[:, i]) for i in 1:size(x, 2)]
    control_interps = [extrapolate(interpolate((s_controls,), u[:, i], Gridded(Linear())), Flat()) for i in 1:size(u, 2)]

    states_interp(t) = [itp(t) for itp in state_interps]
    controls_interp(t) = [itp(t) for itp in control_interps]
    return Result_interpolation(states_interp, controls_interp, s_states)
end

make_result_interpolation(x::AbstractMatrix{Float64}, u::AbstractMatrix{Float64}, s::Vector{Float64}) = make_result_interpolation(x, u, s, s)

const INITIALIZATION_STATE_LABELS = ["vx", "vy", "ψ", "ω", "n", "t"]
const INITIALIZATION_CONTROL_LABELS = ["torque", "steering"]


function build_initialization_controller(car::Car)
    return (
        vref = 5.0,
        steeringP = 10.0,
        steeringD = 1.0,
        velocityP = 3 * car.carParameters.mass.value / 280.0,
    )
end


function get_initialization_controls(car::Car, track::Track, s, x::AbstractVector, controller::NamedTuple)
    max_steer = car.wheelAssemblies[1].maxAngle.value
    th = track.fcurve(s)[2]
    epsilon = x[3] - th
    n_dot = x[1] * sin(epsilon) + x[2] * cos(epsilon)

    torque = (controller.vref - x[1]) * controller.velocityP
    steering = clamp(-(x[5] * controller.steeringP + n_dot * controller.steeringD), -max_steer, max_steer)
    return torque, steering
end


function get_initialization_controls(car::Car, track::Track, s, x::AbstractVector, steeringP, steeringD, velocityP, vref)
    controller = (steeringP=steeringP, steeringD=steeringD, velocityP=velocityP, vref=vref)
    return get_initialization_controls(car, track, s, x, controller)
end


function create_initialization_debug_plot()
    n_state_plots = length(INITIALIZATION_STATE_LABELS)
    n_control_plots = length(INITIALIZATION_CONTROL_LABELS)
    n_plots = n_state_plots + n_control_plots

    debug_fig = Figure(size=(1200, 800))
    debug_axes = [
        Axis(
            debug_fig[div(i - 1, 3) + 1, mod(i - 1, 3) + 1],
            title=(i <= n_state_plots ? INITIALIZATION_STATE_LABELS[i] : INITIALIZATION_CONTROL_LABELS[i - n_state_plots]),
        ) for i in 1:n_plots
    ]
    debug_obs = [Observable(Point2f[]) for _ in 1:n_plots]

    for i in 1:n_plots
        lines!(debug_axes[i], debug_obs[i]; linewidth=2)
    end
    display(GLMakie.Screen(), debug_fig)

    return debug_axes, debug_obs
end


function update_initialization_debug!(debug_axes, debug_obs, s, x::AbstractVector, torque, steering)
    n_state_plots = length(INITIALIZATION_STATE_LABELS)
    for i in 1:n_state_plots
        push!(debug_obs[i][], Point2f(s, x[i]))
        notify(debug_obs[i])
        autolimits!(debug_axes[i])
    end

    control_values = (torque, steering)
    for (j, value) in enumerate(control_values)
        idx = n_state_plots + j
        push!(debug_obs[idx][], Point2f(s, value))
        notify(debug_obs[idx])
        autolimits!(debug_axes[idx])
    end
end


function create_initialization_callback(span2, debug_axes, debug_obs, car::Car, track::Track, controller)
    s_start, s_end = span2
    return FunctionCallingCallback(; funcat=range(s_start, s_end, length=201)) do x, s, integrator
        pct = round(100 * (s - s_start) / (s_end - s_start), digits=0)
        print("\r  Initialization: $(Int(pct))% (s=$(round(s, digits=2))/$(round(s_end, digits=2)))")

        torque, steering = get_initialization_controls(car, track, s, x, controller)
        update_initialization_debug!(debug_axes, debug_obs, s, x, torque, steering)
    end
end


function create_initialization_controls_trajectory(car::Car, track::Track, s::AbstractVector, x::AbstractMatrix, controller)
    torque = zeros(length(s))
    steering = zeros(length(s))
    for i in eachindex(s)
        torque[i], steering[i] = get_initialization_controls(car, track, s[i], view(x, i, :), controller)
    end
    return torque, steering
end


function plot_initialization_summary(track::Track, initialization::Result_interpolation, s::AbstractVector, x::AbstractMatrix)
    fig_init = Figure()
    ax_init = Axis(fig_init[1, 1], aspect=DataAspect(), title="Initialization")
    plotCarPath_interpolated(track, initialization, ax_init)
    display(GLMakie.Screen(), fig_init)

    state_history_labels = ["vx", "vy", "ψ", "ψ̇", "n", "t"]
    fig_states = Figure()
    for j in 1:size(x, 2)
        ax = Axis(fig_states[j, 1], ylabel=state_history_labels[j])
        lines!(ax, s, x[:, j], linewidth=2)
    end
    display(GLMakie.Screen(), fig_states)
end


function unpack_initialization_params(p)
    if p isa NamedTuple
        return p.car, p.track, p.controller
    end

    controller = (steeringP=p[3], steeringD=p[4], velocityP=p[5], vref=p[6])
    return p[1], p[2], controller
end


function initializeSolution_interpolation(car::Car, track::Track, segments::Int64)
    println("started initialization")
    controller = build_initialization_controller(car)

    x0 = [controller.vref, 0.0, track.theta[1], 0.0, 0.0, 0.0]
    sampling_distances = LinRange(track.sampleDistances[1], track.sampleDistances[end], segments)
    span2 = (sampling_distances[1], sampling_distances[end])
    problem_params = (car=car, track=track, controller=controller)
    prob = ODEProblem(carODE_path_initialization, x0, span2, problem_params)

    debug_axes, debug_obs = create_initialization_debug_plot()
    cb = create_initialization_callback(span2, debug_axes, debug_obs, car, track, controller)

    sol = OrdinaryDiffEq.solve(prob, Rodas4(autodiff = AutoFiniteDiff()), saveat=sampling_distances, reltol=1e-1, abstol=1e-1, callback=cb)
    println()

    x = hcat(sol.u...)'
    s = sol.t
    torque, steering = create_initialization_controls_trajectory(car, track, s, x, controller)

    u = zeros(segments, Int64(car.carParameters.nControls.value))
    u[:, 1] = torque
    u[:, 2] = steering

    initialization = make_result_interpolation(x, u, s)
    plot_initialization_summary(track, initialization, s, x)

    return initialization
end;


function carODE_path(car::Car, track::Track, k::Union{Int64,Float64},
    u::AbstractVector, x::AbstractVector,
    model::Union{JuMP.Model,Nothing}=nothing)
    #car.mapping(car,u,x)
    car.controlMapping(u)
    car.stateMapping(x)
    dzds = time2path(car, track, k, model) #time2path(s,instantCarParams,track,car)
    return dzds
end;

#poradie car,track,k,s,u,x


function carODE_path_initialization(du, x, p, s)
    car, track, controller = unpack_initialization_params(p)
    torque, steering = get_initialization_controls(car, track, s, x, controller)
    control = [torque, steering, 0.0]
    dx = carODE_path(car, track, s, control, x, nothing) # carF must return a vector matching length(x)
    du .= dx
end;


#initializeSolution(1,2)
#here I need to define transfomation of ODE with respect to time to ODE with respect to path
function time2path(car::Car, track::Track, k::Union{Int64,Float64}, model::Union{Nothing,JuMP.Model})
    #track.mapping(track,instantTrack,s)
    dxdt = car.carFunction(track, model)
    v_x = car.carParameters.velocity.value[1]
    v_y = car.carParameters.velocity.value[2]
    psi = car.carParameters.psi.value
    n = car.carParameters.n.value

    th = track.fcurve(k)[2]
    C = track.fcurve(k)[1]

    epsilon = psi - th

    Sf = (1 - n * C) ./ (v_x .* cos(epsilon) - v_y .* sin(epsilon))
    dndt = v_x .* sin(epsilon) + v_y .* cos(epsilon)

    dzds = [
        Sf .* dxdt[1];    #dvx/ds
        Sf .* dxdt[2];    #dvy/ds
        Sf .* dxdt[3];    #dpsi/ds
        Sf .* dxdt[4];    #ddpsi/ds
        Sf .* dndt;       #dn/ds
        Sf               #dt/ds 
        # this should be capable of accepting more states
    ]

    #Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
    return dzds
end

function findOptimalTrajectory(track::Track, car::Car, model::JuMP.Model, sampleDistances, initialization=nothing)
    N = length(sampleDistances)

    #fill track parameters which are constant along track, should be automatized
    track.rho = fill(track.rho[1], N)
    track.μ = fill(track.μ[1], N)

    nStates = Int(car.carParameters.nStates.value)
    s = sampleDistances
    #create inputs for car model #create instance of track parameters

    #determine sizes of inputs and states
    nControls = Int(round(car.carParameters.nControls.value))

    function f(x, u, s, model)
        dxds = carODE_path(car, track, s, u, x, model)
        return dxds
    end

    stages = 2
    lobotom = createLobattoIIIA(stages, f)
    xd = lobotom.createConstraints(f, 6, 3, track.fcurve, s, model, initialization.states, initialization.controls)
    X = xd[2]
    U = xd[3][1:end-1, :]
    s_all = xd[6]
    @objective(model, Min, X[end, 6])


    @constraint(model, X[1:end, 1] .>= 0) #vx
    #@constraint(model,X[1,1] .== 5) # intial vx
    @constraint(model, X[1, 2] .== 0) # intial vy
    @constraint(model, X[1, 3] .== track.theta[1]) # intial heading
    #@constraint(model,X[end,3].== track.theta[end])
    #@constraint(model,X[end,5].== 0)
    @constraint(model, X[1, 6] .>= 0) # final time
    @constraint(model, diff(X[:, 6]) .>= 0) #time goes forward



    optimize!(model)

    x = value.(X)
    u = value.(U)
    out = Result(x, u, s_all)
    out_interp = make_result_interpolation(x, u, s_all)
    return out, out_interp
end

function find_optimal_trajectory2(problem::Problem_config, segments::Int64, pol_order::Int64, variant)
    track = problem.track
    car = problem.car
    model = problem.model
    function F(x, u, s)
        #@infiltrate
        dxds = carODE_path(car, track, s, u, x, model)
        return dxds
    end

    function f(x, u, s, model)
        dxds = carODE_path(car, track, s, u, x, model)
        return dxds
    end

    nControls = Int64(car.carParameters.nControls.value)
    nStates = Int64(car.carParameters.nStates.value)
    initialization = initializeSolution_interpolation(car, track, 200)

    Gauss_radau = create_gauss_pseudospectral_metod(F,pol_order,variant,model,nControls,nStates,track);
    (params, tunables) = setParameters(car,model)
    problem.params = params
    xd = Gauss_radau.createConstraints(segments,initialization);

    #Gauss_legendre = create_gauss_legendre(F,pol_order,variant,model,nControls,nStates,track);
    #xd = Gauss_legendre.createConstraints(segments,initialization);
    #segment_edges = LinRange(track.sampleDistances[1], track.sampleDistances[end], segments + 1)
    #RKadaptive = createLobattoIIIA_Adaptive(f, pol_order, model, nControls, nStates, track)
    #xd = RKadaptive.createConstraints(segment_edges, initialization)
    X = xd[2]
    U = xd[3]
    s_all = xd[4]
    segment_edges = xd[5]

    @constraint(model, X[1:end, 1] .>= 0) #vx
    @constraint(model, X[1, 1] .>= 2) #vx
    @constraint(model, X[1, 2] .== 0) # intial vy
    #@constraint(model, X[1, 3] .== track.theta[1]) # intial heading
    @constraint(model, X[1, 6] .>= 0) # final time
    @constraint(model, diff(X[:, 6]) .>= 0) #time goes forward


    #@constraint(model,diff(U[:,1]) .>=-1) #constraint on controls derivative

    #this has to have variable length to allow closed track
    #@constraint(model,-100/2 .<= diff(U[1:end,1])./diff(X[:,6]) .<= 100/2) #constraint on controls derivative
    #@constraint(model,-100/2 .<= diff(U[1:end,2])./diff(X[:,6]) .<= 100/2) #constraint on controls derivative
    #@constraint(model,-1/2 .<= diff(U[1:end-1,3]./diff(X[:,6])) .<= 1/2) #constraint on controls derivative
    #@constraint(model,U[1,:].== U[2,:])
    @objective(model, Min, X[end, 6])
    optimize!(model)
    resetParameters(tunables)
    x = value.(X)
    u = value.(U)

    # Interpolate using the LGR polynomial from the optimisation (Garg et al. 2010)
    #out_interp = Gauss_legendre.createInterpolator(x, u, s_all, segment_edges)
    out_interp = Gauss_radau.createInterpolator(x, u, s_all, segment_edges)
    #out_interp = RKadaptive.createInterpolator(x, u, s_all, segment_edges)

    out = Result(x, u[1:end-1, :], s_all)
    return out, out_interp
end

function refineMesh(problem,segment_edges,s_all,pol_order)
    error_itp = getErrors(problem)
    segment_errors = zeros(length(segment_edges) - 1)
    idx = 1
    clear = 1
    for i = eachindex(segment_errors)
        error_segment = 0
        for node = 2:pol_order
            error_segment += error_itp(s_all[idx])
            idx += 1
        end
        segment_errors[i] = error_segment
    end
    error_threshold = 1e-0

    # Find and insert nodes for segments with error > 1e-3
    segment_edges = collect(segment_edges)  # Convert LinRange to Vector for insertion
    for i = length(segment_errors):-1:1  # Iterate backwards to avoid index shifting issues
        if segment_errors[i] > error_threshold
            s_mid = (segment_edges[i] + segment_edges[i+1]) / 2
            insert!(segment_edges, i + 1, s_mid)
            clear = 0
        end
    end
    fig = plot(segment_errors, axis=(; yscale=log10))
    lines!(fig.axis, [1, length(segment_errors)], [error_threshold, error_threshold], color=:red, linewidth=2)
    display(fig)
    sleep(1)   # give the backend time to render
    return segment_edges,clear,segment_errors
end




function find_optimal_trajectory_adaptive(problem::Problem_config, segments::Int64, pol_order::Int64, variant)
    track = problem.track
    car = problem.car
    model = problem.model

    nControls = Int64(car.carParameters.nControls.value)
    nStates = Int64(car.carParameters.nStates.value)
    x = nothing
    u = nothing
    s_all = nothing

    function f(x, u, s, model)
        dxds = carODE_path(car, track, s, u, x, model)
        return dxds
    end

    initialization = initializeSolution_interpolation(car, track, 200)
    segment_edges = LinRange(track.sampleDistances[1], track.sampleDistances[end], segments + 1)
    clear = 0
    iterations = 0
    while clear == 0 && iterations <= 50
        clear = 1
        iterations += 1
        RKadaptive = createLobattoIIIA_Adaptive(f, pol_order, model, nControls, nStates, track)
        println("creating constraints")

        #Add parameters
        (params, tunables) = setParameters(car,model)
        problem.params = params
        xd = RKadaptive.createConstraints(segment_edges, initialization)
        X = xd[2]
        U = xd[3]
        s_all = xd[4]
        segment_edges = xd[5]

        @constraint(model, X[1:end, 1] .>= 0) #vx
        @constraint(model, X[1, 1] .== 10) #vx
        @constraint(model, X[1, 2] .== 0) # intial vy
        @constraint(model, X[1, 3] .== track.theta[1]) # intial heading
        @constraint(model, X[1, 6] .>= 0) # final time
        @constraint(model, diff(X[:, 6]) .>= 0) #time goes forward
        @objective(model, Min, X[end, 6])

        #@constraint(model,-50 .<= diff(U[1:end,1])./diff(X[:,6]) .<= 50) #constraint on controls derivative
        #@constraint(model,-50 .<= diff(U[1:end,2])./diff(X[:,6]) .<= 50) #constraint on controls derivative
        #@constraint(model,-0.5 .<= diff(U[1:end-1,3]./diff(X[:,6])) .<= 0.5) #constraint on controls derivative

        println("Variables: $(num_variables(model)), Constraints: $(num_constraints(model; count_variable_in_set_constraints=true))")
        set_optimizer_attribute(model, "max_iter", 3000)
        set_optimizer_attribute(model, "print_level", 5)
        set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
        set_optimizer_attribute(model, "mu_strategy", "adaptive")
        optimize!(model)
        println("Termination: ", termination_status(model), " | Objective: ", objective_value(model))

        x = value.(X)
        u = value.(U)
        problem.model = model
        problem.params = params
        problem.optiResult = RKadaptive.createInterpolator(x, u, s_all, segment_edges)
        #initialization = make_result_interpolation(x, u, s_all) # uncomment to use previous result as initialization

        resetParameters(tunables)
        (segment_edges,clear,segment_errors) = refineMesh(problem,segment_edges,s_all,pol_order)

        if clear == 0
            model = DiffOpt.nonlinear_diff_model(Ipopt.Optimizer)
        end

        println("Sum of all errors: ", sum(segment_errors))


    end
    out = Result(x, u[1:end-1, :], s_all)
    return out, problem.optiResult
end