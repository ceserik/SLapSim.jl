
function refineMesh(problem, segment_edges, pol_order; iteration::Int=0)
    cfg = problem.mesh_refinement
    segment_errors = getSegmentErrors(problem; method=cfg.error_method)
    clear = 1
    any_bad = false

    segment_edges = collect(segment_edges)
    for i = length(segment_errors):-1:1
        if segment_errors[i] > cfg.tol
            any_bad = true
            if cfg.method === :h || cfg.method === :hp
                s_mid = (segment_edges[i] + segment_edges[i+1]) / 2
                insert!(segment_edges, i + 1, s_mid)
            end
            clear = 0
        end
    end
    if any_bad && (cfg.method === :p || cfg.method === :hp)
        pol_order += 1
        println("p-refinement: pol_order → $pol_order")
    end
    fig = plot(segment_errors, axis=(; yscale=log10, title="Segment errors — iter $iteration (p=$pol_order)"))
    lines!(fig.axis, [1, length(segment_errors)], [cfg.tol, cfg.tol], color=:red, linewidth=2)
    display(GLMakie.Screen(), fig)
    
    return segment_edges, clear, segment_errors, pol_order
end

struct Result
    states::Matrix{Float64}
    controls::Matrix{Float64}
    path::Vector{Float64}
end


struct Result_interpolation
    states
    controls
    path
    segment_edges
    state_deriv      # s -> dx/ds via analytic differentiation of the collocation polynomial; nothing if unavailable
end
Result_interpolation(states, controls, path, segment_edges) =
    Result_interpolation(states, controls, path, segment_edges, nothing)

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


function find_optimal_trajectory_adaptive(exp::Experiment, segments::Int64, pol_order::Int64, variant)
    track = exp.track
    car = exp.car
    model = exp.model
    performSensitivity = wants_sensitivity(exp.solver)

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
    initialization = initializeSolution_interpolation(car, track, Int64(round(track.sampleDistances[end] * 2)); plot_initialization=exp.analysis.plot_initialization)
    # Create sampling points for track
    segment_edges = LinRange(track.sampleDistances[1], track.sampleDistances[end], segments + 1)

    #From track that has data in vectors: vector of X, vector of Y vector of curvature ... create track that uses functions to describe its states
    # the interpoation function should be from the smoothing, but i am using basic interpolation of julia, work for future
    track_interp = interpolate_track(track)

    #this is used  to create constraints on state n, which gets its constraints from track width
    n_state_index = findfirst(e -> e.name == "n", car.carParameters.state_descriptor)

    clear = 0
    iterations = 0
    while clear == 0 && iterations <= exp.mesh_refinement.max_iterations
        clear = 1
        iterations += 1
        #@infiltrate

        # Scaling of controls and states, scaling is derived from bounds setup when creating a car
        # These functions 
        x_scale = get_scales(car.carParameters.state_descriptor)
        u_scale = get_scales(car.carParameters.control_descriptor)
        x_lb_static, x_ub_static = get_bounds(car.carParameters.state_descriptor)
        u_lb_static, u_ub_static = get_bounds(car.carParameters.control_descriptor)

        margin = car.chassis.track.value/2
        
            n_scale = max(maximum(track.widthL), maximum(track.widthR))
            x_scale[n_state_index] = n_scale
            x_lb_static[n_state_index] = -maximum(track.widthR) + margin
            x_ub_static[n_state_index] = maximum(track.widthL) - margin
        
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

        sp = scaledProblem(f, x_lb_fun, x_ub_fun, u_lb_fun, u_ub_fun, initialization, x_scale, u_scale)
        printBounds(sp)
        RKadaptive = if variant == "Legendre"
            create_gauss_legendre(sp.f, pol_order, variant, model, nControls, nStates, track)
        elseif variant == "Radau"
            create_radau(sp.f, pol_order, model, nControls, nStates, track)
        elseif variant == "LobattoIIIA_integral"
            createLobattoIIIA_Adaptive(sp.f, pol_order, model, nControls, nStates, track)
        else
            error("Unknown variant: $variant. Use \"Legendre\", \"Radau\", or \"LobattoIIIA_integral\".")
        end
        println("creating constraints")

        # Add parameters for sensitivity analysis, makes the problem little bigger, because parameters of car are now variables that ipopt sees
        sensitivityParams = carParameter[]
        if performSensitivity
            (params, sensitivityParams) = setParameters(car, model)
            exp.params = params
        end

        # Promote :design carParameters to free scalar JuMP variables (one per param, constant across horizon)
        (designRefs, designParams) = setDesignVariables(car, model)
        exp.designRefs = designRefs
        xd = RKadaptive.createConstraints(segment_edges, sp.init)
        X_s = xd[2]
        U_s = xd[3]
        s_all = xd[4]
        segment_edges = xd[5]
        s_controls = xd[6]
        applyBounds!(sp, X_s, U_s, s_all, s_controls)
        #scale back
        X = X_s .* sp.x_scale'
        U = U_s .* sp.u_scale'

        # Boundary conditions are supplied by the discipline.
        apply_boundary_conditions!(exp.discipline, model, X, track)

        # Global (lap-integral) constraints from the experiment.
        for gc in exp.global_constraints
            apply_global!(gc, model, X, U, car, s_all)
        end

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
        exp.model = model
        # The resulting control is not just a vector for each state, but it is a function of path with fancy interpolation depending on discretization method
        exp.optiResult = RKadaptive.createInterpolator(x, u, s_all, segment_edges)

        # Here you can use the result of previous optimisiation as warm start, sometimes works very good sometimes very bad
        #initialization = make_result_interpolation(x, u, s_all)
        if performSensitivity
            println("resetting parameters for sensitivity analysis")
            exp.params = params
            resetParameters(sensitivityParams)
        end

        # Write optimized design var values back into carParameter.value (replace JuMP refs with floats)
        if !isempty(designParams)
            resolveDesignVariables(designParams)
            for p in designParams
                println("design var $(p.name) = $(p.value) [$(p.unit)]")
            end
        end

        # Call to Mesh refinment algorigth, if clear == 1 the solution does not need to be refined and optimization ends
        (segment_edges, clear, segment_errors, pol_order) = refineMesh(exp, segment_edges, pol_order; iteration=iterations)

        if clear == 0
            #create clear model for next iteration
            model = build_model(exp.solver)
        end
        println("Sum of all errors: ", sum(segment_errors))
    end
    out = Result(x, u[1:end-1, :], s_all)
    return out, exp.optiResult
end