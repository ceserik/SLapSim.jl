using Revise
using Infiltrator
using SLapSim
using JuMP


struct Initialization
    states::Matrix{Float64}
    controls::Matrix{Float64}
    path::Vector{Float64}
end


function carODE_path(car::Car,track::Track,k::Union{Int64,Float64},u::Union{Vector{VariableRef},Vector{Float64}},x::Union{Vector{VariableRef},Vector{Float64}},model::Union{JuMP.Model,Nothing})
    #car.mapping(car,u,x)
    car.controlMapping(car,u)
    car.stateMapping(car,x)
    #@infiltrate
    dzds = time2path(car,track,k,model) #time2path(s,instantCarParams,track,car)
    return dzds
end;

#poradie car,track,k,s,u,x



function initializeSolution(car::Car,track::Track)
    span2 = [track.sampleDistances[1],track.sampleDistances[end]]
    x0 = [2.0, 0.0, track.theta[1], 0.0, 0, 0.0]
    steeringP = 10
    velocityP = 3
    vref = 2

    prob = ODEProblem(carODEPath, x0, span2,[car,track,steeringP,velocityP,vref]);
    sol = solve(prob, Tsit5(),saveat=track.sampleDistances);

    #@infiltrate
    x = hcat(sol.u...)'
    s = sol.t
    steering = x[:,5] .* steeringP
    torque = (vref .- x[:,1]) * velocityP

    u = zeros(length(track.sampleDistances),Int64(car.carParameters.nControls.value))
    u[:,1] = torque
    u[:,3] = steering
    #@infiltrate
    initialization = Initialization(
        x,
        u[1:end-1,:],
        s
    )
    return initialization
end;


function carODEPath(du, x, p, s)
    car = p[1];
    track = p[2]
    steeringP = p[3]
    velocityP = p[4]
    vref = p[5]
    
    control = [(vref - car.carParameters.velocity.value[1])*velocityP, 0.0, x[5]*steeringP]

    dx = carODE_path(car,track,s, control, x,nothing)  ; # carF must return a vector matching length(x)
    du .= dx;
end;


#initializeSolution(1,2)
## here I need to define trnasofmation of ODE with respect to time to ODE with respect to path
function time2path(car::Car,track::Track,k::Union{Int64,Float64},model::Union{Nothing,JuMP.Model})
    #track.mapping(track,instantTrack,s)
    dxdt = car.carFunction(car,track,k,model)
    v_x = car.carParameters.velocity.value[1]
    v_y = car.carParameters.velocity.value[2]
    psi = car.carParameters.psi.value
    n   = car.carParameters.n.value
    
    if typeof(k) == Float64
        th = interp1(track.sampleDistances, track.theta, k)
        C  = interp1(track.sampleDistances, track.curvature, k)
    else
        th = track.theta[k]
        C = track.curvature[k]
    end
    # until better track model is added epsilon = 0
    #epsilon = psi-th #where psi is heading of car realtive to intertial frame, th is heading of track

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
    #@infiltrate
    #Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
    return dzds
end




function findOptimalTrajectory(track::Track,car::Car,model::JuMP.Model, initialization=nothing)
    N = length(track.sampleDistances)

    #fill track parameters which are constant along track, should be automatized
    track.rho    = fill(track.rho[1],N)
    track.μ      = fill(track.μ[1],N)
    track.widthL = fill(track.widthL[1],N)
    track.widthR = fill(track.widthR[1],N)

    nStates =  Int(car.carParameters.nStates.value)
    #create inputs for car model #create instance of track parameters

    #determine sizes of inputs and states
    nControls = Int(round(car.carParameters.nControls.value))
    @variable(model,U[1:N-1,1:nControls])

    #stateInit = [fill(5,N),zeros(N),zeros(N),zeros(N),zeros(N),range(0,4,N)]'
    @variable(model,X[1:N,1:6])


    if isnothing(initialization)
        for i = 1:N
            set_start_value(X[i,1], 5.0)  # vx = 5
            set_start_value(X[i,2], 0.0)  # vy = 0
            set_start_value(X[i,3], pi/2)  # psi = 0
            set_start_value(X[i,4], 0.0)  # dpsi = 0
            set_start_value(X[i,5], 0.0)  # n = 0
            set_start_value(X[i,6], 4*(i-1)/(N-1))  # t = linspace(0,4,N)
            #println(4*(i-1)/(N-1))
        end
    else
        x = initialization.states
        u = initialization.controls
        s = initialization.path
        for i = 1:N
            set_start_value(X[i,1], x[i,1])  # vx from initialization
            set_start_value(X[i,2], x[i,2])  # vy from initialization
            set_start_value(X[i,3], x[i,3])  # psi from initialization
            set_start_value(X[i,4], x[i,4])  # dpsi from initialization
            set_start_value(X[i,5], x[i,5])  # n from initialization
            set_start_value(X[i,6], x[i,6])  # t from initialization
        end
        for i = 1:N-1
            for j = 1:nControls
                set_start_value(U[i,j], u[i,j])  # controls from initialization
            end
        end
    end
    
    s = track.sampleDistances#1:length(track.curvature)
    #t_final = X[end,6]


    for k = 1:N-1 # loop over control intervals
        x_next = X[k,:] .+ (s[k+1]-s[k]) .* carODE_path(car,track,k, U[k,:], X[k,:],model);
        @constraint(model,X[k+1,:] .== x_next)
    end

    #initial states
    #        [vx vy psi dpsi n t]
    #initial = [4; 0; 0; 0; 0; 0]
    #set initial states
    #@constraint(model,X[1,:] .== initial)

    @constraint(model,X[1:end,1] .>= 0) #vx
    @constraint(model,X[1,1] .== 5) # intial vx
    @constraint(model,X[1,2] .== 0) # intial vx
    @constraint(model,X[1,6] .>= 0) # final time
    @constraint(model,diff(X[:,6]) .>=0) #time goes forward


    @objective(model,Min,X[end,6])
    optimize!(model)

    # Plotting the results
    fig = Figure(layout = GridLayout(6, 1))
    
    # Create subplots for each variable
    labels = ["vx", "vy", "psi", "dpsi", "n", "t"]
    for (i, label) in enumerate(labels)
        ax = Axis(fig[i, 1], ylabel = label)
        lines!(ax, value.(X[:,i]), label=label)
    end
    display(GLMakie.Screen(), fig)  # This creates a new window

    out = Initialization(value.(X),value.(U),track.sampleDistances)
    return out
end

