
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

function initializeSolution_interpolation(car::Car, track::Track, segments::Int64; vref=5.0)
    println("started initialization")
    max_steer = car.wheelAssemblies[1].maxAngle.value
    Kv = 3 * car.carParameters.mass.value / 280.0
    Kp, Kd = 8.0, 3.0
    nControls = Int64(car.carParameters.nControls.value)

    function ctrl(s, x)
        th = track.fcurve(s)[2]
        ε = x[3] - th
        n_dot = x[1] * sin(ε) + x[2] * cos(ε)
        torque = (vref - x[1]) * Kv
        steering = clamp(-(x[5] * Kp + n_dot * Kd), -max_steer, max_steer)
        return torque, steering
    end

    x0 = [vref, 0.0, track.theta[1], 0.0, 0.0, 0.0]
    s_span = (track.sampleDistances[1], track.sampleDistances[end])
    s_save = LinRange(s_span..., segments)

    # Live debug plot
    labels = ["vx", "vy", "ψ", "ψ̇", "n", "t", "torque", "steering"]
    fig = Figure(size=(1200, 800))
    debug_axes = [Axis(fig[div(i-1, 3)+1, mod(i-1, 3)+1], title=labels[i]) for i in eachindex(labels)]
    debug_obs = [Observable(Point2f[]) for _ in eachindex(labels)]
    for i in eachindex(labels)
        lines!(debug_axes[i], debug_obs[i]; linewidth=2)
    end
    display(GLMakie.Screen(), fig)

    cb = FunctionCallingCallback(; funcat=range(s_span..., length=201)) do x, s, _
        pct = round(100 * (s - s_span[1]) / (s_span[2] - s_span[1]), digits=0)
        print("\r  Initialization: $(Int(pct))%")
        torque, steering = ctrl(s, x)
        for i in 1:6
            push!(debug_obs[i][], Point2f(s, x[i])); notify(debug_obs[i]); autolimits!(debug_axes[i])
        end
        for (j, v) in enumerate((torque, steering))
            push!(debug_obs[6+j][], Point2f(s, v)); notify(debug_obs[6+j]); autolimits!(debug_axes[6+j])
        end
    end

    prob = ODEProblem((du, x, _, s) -> begin
        torque, steering = ctrl(s, x)
        du .= carODE_path(car, track, s, [torque, steering, zeros(nControls - 2)...], x, nothing)
    end, x0, s_span)

    sol = OrdinaryDiffEq.solve(prob, Rodas4(autodiff=AutoFiniteDiff()), saveat=s_save, reltol=1e-4, abstol=1e-4, callback=cb)
    println()

    x = hcat(sol.u...)'
    s = sol.t
    u = zeros(segments, nControls)
    for i in eachindex(s)
        u[i, 1], u[i, 2] = ctrl(s[i], view(x, i, :))
    end

    initialization = make_result_interpolation(x, u, s)

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
    control_reg = 1e-6 * sum(U[i,j]^2 for i in axes(U,1), j in axes(U,2))
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
    initialization = initializeSolution_interpolation(car, track, 200)

    Gauss_legendre = create_gauss_legendre(F,pol_order,variant,model,nControls,nStates,track);
    #Gauss_radau = create_gauss_pseudospectral_metod(F,pol_order,variant,model,nControls,nStates,track);
    (params, tunables) = setParameters(car,model)
    problem.params = params
    
    
    #xd = Gauss_radau.createConstraints(segments,initialization);
    xd = Gauss_legendre.createConstraints(segments,initialization);
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


    #@constraint(model,diff(U[:,1]) .>=-1) #constraint on controls derivative

    #this has to have variable length to allow closed track
    #@constraint(model,-100/2 .<= diff(U[1:end,1])./diff(X[:,6]) .<= 100/2) #constraint on controls derivative
    #@constraint(model,-100/2 .<= diff(U[1:end,2])./diff(X[:,6]) .<= 100/2) #constraint on controls derivative
    #@constraint(model,-1/2 .<= diff(U[1:end-1,3]./diff(X[:,6])) .<= 1/2) #constraint on controls derivative
    #@constraint(model,U[1,:].== U[2,:])
    control_reg = 1e-6 * sum(U[i,j]^2 for i in axes(U,1), j in axes(U,2))
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
    error_threshold = 1e-1

    # Find and insert nodes for segments with error > 1e-3
    segment_edges = collect(segment_edges)  # Convert LinRange to Vector for insertion
    for i = length(segment_errors):-1:1  # Iterate backwards to avoid index shifting issues
        if segment_errors[i] > error_threshold
            s_mid = (segment_edges[i] + segment_edges[i+1]) / 2
            insert!(segment_edges, i + 1, s_mid)
            clear = 0
        end
    end
    fig = plot(segment_errors, axis=(; yscale=log10, title="Segment errors"))
    lines!(fig.axis, [1, length(segment_errors)], [error_threshold, error_threshold], color=:red, linewidth=2)
    display(GLMakie.Screen(), fig)
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
        # Characteristic scales so optimizer variables are O(1)
        x_scale = [10.0, 5.0, 1.0, 1.0, 5.0, 10.0]  # [vx, vy, ψ, ω, n, t]
        u_scale = [30.0, 0.3, 30.0]                    # [torque, steering, torque]
        RKadaptive = createLobattoIIIA_Adaptive(f, pol_order, model, nControls, nStates, track;
                                                 x_scale=x_scale, u_scale=u_scale)
        println("creating constraints")

        #Add parameters
        (params, tunables) = setParameters(car,model)
        problem.params = params
        xd = RKadaptive.createConstraints(segment_edges, initialization)
        X = xd[2]
        U = xd[3]
        s_all = xd[4]
        segment_edges = xd[5]

        # State bounds: [vx, vy, ψ, ω, n, t]
        @constraint(model, X[1, 1] .== 10)      # vx
        @constraint(model, X[1, 2] .>= 0)    # vy
        @constraint(model, X[1, 3] .== track.theta[1]) # intial heading
        @constraint(model, X[1, 4] .== 0)     # ω
        @constraint(model, X[1, 6] .== 0)      # t
        @constraint(model, X[1, 2] .== 0) # intial vy
        
        @constraint(model, X[1, 6] .>= 0) # final time
        @constraint(model, diff(X[:, 6]) .>= 0) #time goes forward
        control_reg = 1e-6 * sum((U[i,j]/u_scale[j])^2 for i in axes(U,1), j in axes(U,2))
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