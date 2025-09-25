using SLapSim
using JuMP
using Infiltrator
function carODE_path(car::Car,track::Track,k::Int64,u::Union{Vector{VariableRef},Vector{Float64}},x::Union{Vector{VariableRef},Vector{Float64}},model::Union{JuMP.Model,Nothing})
    #car.mapping(car,u,x)
    car.controlMapping(car,u)
    car.stateMapping(car,x)
    dzds = time2path(car,track,k,model) #time2path(s,instantCarParams,track,car)
    return dzds
end;

#poradie car,track,k,s,u,x



function initializeSolution(car::Car,track::Track)
    span = [track.sampleDistances[1],track.sampleDistances[end]]

    prob = ODEProblem(f!, x0, tspan,car);
    sol = solve(prob, Tsit5(), saveat=0.01);
    return span
end;


function f!(du, x, p, t)
    car = p[1];
    track = p[2]
    k = p[3]
    control = u_const(t);
    dx = carODE_path(car,track,k, control, x,nothing)  ; # carF must return a vector matching length(x)
    #print(dx, "\n")
    du .= dx;
end;


#initializeSolution(1,2)
## here I need to define trnasofmation of ODE with respect to time to ODE with respect to path
function time2path(car::Car,track::Track,k::Int64,model::Union{Nothing,JuMP.Model})
    #track.mapping(track,instantTrack,s)
    dxdt = car.carFunction(car,track,k,model)
    v_x = car.carParameters.velocity.value[1]
    v_y = car.carParameters.velocity.value[2]
    psi = car.carParameters.psi.value
    n   = car.carParameters.n.value

    th = track.theta[k]
    C = track.curvature[k]
    # until better track model is added epsilon = 0
    #epsilon = psi-th #where psi is heading of car realtive to intertial frame, th is heading of track

    epsilon = psi - th

    Sf = (1 - n*C) ./ (v_x.*cos(epsilon) - v_y.*sin(epsilon));
    dndt = v_x.*sin(epsilon) + v_y.*cos(epsilon);

    dzds = [
    Sf .* dxdt[1];    #dvxds
    Sf .* dxdt[2];    #dvyds
    Sf .* dxdt[3];    #dpsids
    Sf .* dxdt[4];    #ddpsids
    Sf .* dndt;       #dnds
    Sf;               #ds
    # this should be capable of accepting more states
    ]
    return dzds
end




function findOptimalTrajectory(track::Track,car::Car,model::JuMP.Model)
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

    for i = 1:N
        set_start_value(X[i,1], 5.0)  # vx = 5
        set_start_value(X[i,2], 0.0)  # vy = 0
        set_start_value(X[i,3], pi/2)  # psi = 0
        set_start_value(X[i,4], 0.0)  # dpsi = 0
        set_start_value(X[i,5], 0.0)  # n = 0
        set_start_value(X[i,6], 4*(i-1)/(N-1))  # t = linspace(0,4,N)
        #println(4*(i-1)/(N-1))
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
    return X, U
end

