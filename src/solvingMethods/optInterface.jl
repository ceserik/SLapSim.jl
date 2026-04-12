
using SLapSim
using JuMP
using OrdinaryDiffEq
using DifferentiationInterface
using Interpolations
using HSL_jll

function make_ipopt_model()
    model = DiffOpt.nonlinear_diff_model(Ipopt.Optimizer)
    #set_attribute(model, "hsllib", HSL_jll.libhsl_path)
    #set_attribute(model, "linear_solver", "ma97")
    #model = Model(Ipopt.Optimizer)
    #set_attribute(model, "pardisolib", "/home/riso/Downloads/panua-pardiso-20240229-linux/lib/libpardiso.so")
    #set_attribute(model, "linear_solver", "pardiso")


    #JuMP.set_optimizer_attribute(model, "max_iter", 3000)
    JuMP.set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
    #JuMP.set_optimizer_attribute(model, "print_timing_statistics", "yes")
    for (k, v) in [
        ("ma57_pre_alloc",     100.0),   
        ("mu_strategy", "adaptive") ,
        ("alpha_for_y",                      "safer-min-dual-infeas"),
        ("recalc_y",                         "yes"),
        ("recalc_y_feas_tol",                1e-4),
        ("adaptive_mu_globalization",        "kkt-error"),
        ("quality_function_balancing_term",  "cubic"),
        ("mu_min",                           1e-10),
        ("nlp_scaling_constr_target_gradient", 1.0),
        ("nlp_scaling_obj_target_gradient",  1.0),
        ("nlp_scaling_min_value",            1e-6),
        ("jacobian_regularization_value",    1e-6),
        ("bound_relax_factor",               1e-6),
        ("mumps_pivtol",                     1e-4),
        ("min_refinement_steps",             4),
        ("max_refinement_steps",             100),
        ("acceptable_dual_inf_tol",          1e-1),
    ]
        JuMP.set_optimizer_attribute(model, k, v)
    end
    return model
end

function refineMesh(problem, segment_edges, s_all, pol_order; error_method::Symbol=:ode)
    error_itp = getErrors(problem; method=error_method)
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
    error_threshold = 1e-1

    # Find and insert nodes for segments with error > 1e-3
    segment_edges = collect(segment_edges) 
    for i = length(segment_errors):-1:1  
        if segment_errors[i] > error_threshold
            s_mid = (segment_edges[i] + segment_edges[i+1]) / 2
            insert!(segment_edges, i + 1, s_mid)
            clear = 0
        end
    end
    fig = plot(segment_errors, axis=(; yscale=log10, title="Segment errors"))
    lines!(fig.axis, [1, length(segment_errors)], [error_threshold, error_threshold], color=:red, linewidth=2)
    display(GLMakie.Screen(), fig)
    return segment_edges, clear, segment_errors
end

struct Result
    states::Matrix{Float64}
    controls::Matrix{Float64}
    path::Vector{Float64}
end


mutable struct Problem_config
    car::Union{Nothing,Car}
    track::Union{Nothing,Track}
    model::Union{Nothing,JuMP.Model}
    optiResult
    params
end

struct Result_interpolation
    states
    controls
    path
    segment_edges
end

function make_result_interpolation(x::AbstractMatrix{Float64}, u::AbstractMatrix{Float64}, s_states::Vector{Float64}, s_controls::Vector{Float64})
    state_interps = [linear_interpolation(s_states, x[:, i]) for i in 1:size(x, 2)]
    control_interps = [extrapolate(interpolate((s_controls,), u[:, i], Gridded(Linear())), Flat()) for i in 1:size(u, 2)]

    states_interp(t) = [itp(t) for itp in state_interps]
    controls_interp(t) = [itp(t) for itp in control_interps]
    return Result_interpolation(states_interp, controls_interp, s_states, Float64[])
end

make_result_interpolation(x::AbstractMatrix{Float64}, u::AbstractMatrix{Float64}, s::Vector{Float64}) = make_result_interpolation(x, u, s, s)


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
    s = sampleDistances
    #determine sizes of inputs and states
    
    nStates = Int(car.carParameters.nStates.value)
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
    control_reg = 1e-6 * sum(U[i, j]^2 for i in axes(U, 1), j in axes(U, 2))
    @objective(model, Min, X[end, 6] + control_reg)


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
    initialization = initializeSolution_interpolation(car, track, 1000)

    Gauss_legendre = create_gauss_legendre(F, pol_order, variant, model, nControls, nStates, track; car=car)
    #Gauss_radau = create_gauss_pseudospectral_metod(F,pol_order,variant,model,nControls,nStates,track);
    (params, tunables) = setParameters(car, model)
    problem.params = params


    #xd = Gauss_radau.createConstraints(segments,initialization);
    xd = Gauss_legendre.createConstraints(segments, initialization)
    X = xd[2]
    U = xd[3]
    s_all = xd[4]
    segment_edges = xd[5]

    @constraint(model, X[1, 1] .<= 10)      # vx
    @constraint(model, X[1, 3] .== track.theta[1]) # intial heading
    @constraint(model, X[1, 4] .== 0)     # ω
    @constraint(model, X[1, 6] .== 0)      # t
    @constraint(model, X[1, 2] .== 0) # intial vy

    @constraint(model, diff(X[:, 6]) .>= 0) #time goes forward
    u_scale = !isnothing(car.control_desc) ? get_scales(car.control_desc) : ones(nControls)
    control_reg = 1e-5 * sum((U[i, j] / u_scale[j])^2 for i in axes(U, 1), j in axes(U, 2))
    @objective(model, Min, X[end, 6] + control_reg)


    optimize!(model)
    resetParameters(tunables)
    x = value.(X)
    u = value.(U)

    # Interpolate using the LGR polynomial from the optimisation (Garg et al. 2010)
    out_interp = Gauss_legendre.createInterpolator(x, u, s_all, segment_edges)
    #out_interp = Gauss_radau.createInterpolator(x, u, s_all, segment_edges)
    #out_interp = RKadaptive.createInterpolator(x, u, s_all, segment_edges)

    out = Result(x, u[1:end-1, :], s_all)
    return out, out_interp
end



function find_optimal_trajectory_adaptive(problem::Problem_config, segments::Int64, pol_order::Int64, variant)
    track = problem.track
    car = problem.car
    model = problem.model

    #get size of control and state vector
    nControls = Int64(car.carParameters.nControls.value)
    nStates = Int64(car.carParameters.nStates.value)
    x = nothing
    u = nothing
    s_all = nothing

    function f(x, u, s, model)
        dxds = carODE_path(car, track, s, u, x, model)
        return dxds
    end

    #initialize solution = drive PD controller along centerline, one controller for all cars right now, but controller should be put inside car structure in future
    initialization = initializeSolution_interpolation(car, track, Int64(round(track.sampleDistances[end] * 2)))
    # Create sampling points for track
    segment_edges = LinRange(track.sampleDistances[1], track.sampleDistances[end], segments + 1)

    #From track that has data in vectors: vector of X, vector of Y vector of curvature ... create track that uses functions to describe its states
    # the interpoation function should be from the smoothing, but i am using basic interpolation of julia, work for future
    track_interp = interpolate_track(track)

    #this is used  to create constraints on state n, which gets its constraints from track width
    n_state_index = findfirst(e -> e.name == "n", car.carParameters.state_descriptor)

    clear = 0
    iterations = 0
    while clear == 0 && iterations <= 50
        clear = 1
        iterations += 1
        #@infiltrate

        # Scaling of controls and states, scaling is derived from bounds setup when creating a car
        # These functions 
        x_scale = get_scales(car.carParameters.state_descriptor)
        u_scale = get_scales(car.carParameters.control_descriptor)
        x_lb_static, x_ub_static = get_bounds(car.carParameters.state_descriptor)
        u_lb_static, u_ub_static = get_bounds(car.carParameters.control_descriptor)

        # state "n" needs special constraints
        margin = car.chassis.track.value/2
        function x_lb_fun(s)
            v = copy(x_lb_static)
            if n_state_index !== nothing
                v[n_state_index] = -track_interp.widthR(s) + margin
            end
            return v
        end
        function x_ub_fun(s)
            v = copy(x_ub_static)
            if n_state_index !== nothing
                v[n_state_index] = track_interp.widthL(s) - margin
            end
            return v
        end
        u_lb_fun(s) = u_lb_static
        u_ub_fun(s) = u_ub_static

        RKadaptive = createLobattoIIIA_Adaptive(f, pol_order, model, nControls, nStates, track;
            x_scale=x_scale, u_scale=u_scale,
            x_lb=x_lb_fun, x_ub=x_ub_fun, u_lb=u_lb_fun, u_ub=u_ub_fun)
        println("creating constraints")

        # Add parameters for sensitivity analysis, makes the problem little bigger, because parameters of car are now variables that ipopt sees
        (params, tunables) = setParameters(car, model)
        problem.params = params
        xd = RKadaptive.createConstraints(segment_edges, initialization)
        X = xd[2]
        U = xd[3]
        s_all = xd[4]
        segment_edges = xd[5]

        # State bounds: [vx, vy, ψ, ω, n, t]
        #these should be depending on scenario: autox, endurance,... I  should create new struct to specify this
        @constraint(model, X[1, 1] .<= 10)      # vx
        @constraint(model, X[1, 3] .== track.theta[1]) # intial heading
        @constraint(model, X[1, 4] .== 0)     # ω
        @constraint(model, X[1, 6] .== 0)      # t
        @constraint(model, X[1, 2] .== 0) # intial vy

        @constraint(model, diff(X[:, 6]) .>= 0) #time goes forward
        control_reg = 0#1e-4 * sum((U[i, j] / u_scale[j])^2 for i in axes(U, 1), j in axes(U, 2))
        @objective(model, Min, X[end, 6] + control_reg)

        #@constraint(model,-50 .<= diff(U[1:end,1])./diff(X[:,6]) .<= 50) #constraint on controls derivative
        #@constraint(model,-50 .<= diff(U[1:end,2])./diff(X[:,6]) .<= 50) #constraint on controls derivative
        #@constraint(model,-0.5 .<= diff(U[1:end-1,3]./diff(X[:,6])) .<= 0.5) #constraint on controls derivative

        println("Variables: $(num_variables(model)), Constraints: $(num_constraints(model; count_variable_in_set_constraints=true))")

        optimize!(model)
        println("Termination: ", termination_status(model), " | Objective: ", objective_value(model))

        x = value.(X)
        u = value.(U)
        problem.model = model
        # The resulting control is not just a vector for each state, but it is a function of path with fancy interpolation depending on discretization method
        problem.optiResult = RKadaptive.createInterpolator(x, u, s_all, segment_edges)

        # Here you can use the result of previous optimisiation as warm start, sometimes works very good sometimes very bad
        #initialization = make_result_interpolation(x, u, s_all)
        println("resetting parameters for sensitivity analysis")
        problem.params = params
        resetParameters(tunables)

        # Call to Mesh refinment algorigth, if clear == 1 the solution does not need to be refind and optimization ends
        (segment_edges, clear, segment_errors) = refineMesh(problem, segment_edges, s_all, pol_order)

        if clear == 0
            #create clear model for next iteration
            model = make_ipopt_model()
        end
        println("Sum of all errors: ", sum(segment_errors))
    end
    out = Result(x, u[1:end-1, :], s_all)
    return out, problem.optiResult
end