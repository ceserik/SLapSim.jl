using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
include("carCreate.jl")
include("trackDefinition.jl")

function f(s,u,x)
    car.mapping(instantCarParams,track,instantTrack,u,x,s)
    dzds = time2path(s,instantCarParams,instantTrack,car)
    return dzds
end

## here I need to define trnasofmation of ODE with respect to time to ODE with respect to path
function time2path(s,instantCarParams,instantTrack,car)
    track.mapping(track,instantTrack,s)
    dxdt = car.carFunction(car,instantCarParams,instantTrack,model)
    v_x = instantCarParams.vx.value
    v_y = instantCarParams.vy.value
    psi = instantCarParams.psi.value
    n = instantCarParams.n.value

    th = track.theta[s]
    C = track.curvature[s]
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

car = createCTU25()
track = singleTurn()

#the track is assumed to be discretized before

model = Model(Ipopt.Optimizer) 
#create inputs for car model #create instance of track parameters
instantCarParams = deepcopy(car.carParameters)
instantTrack = deepcopy(track)

cas = 1000
tfinal_0 = cas
#determine sizes of inputs and states
N = length(track.curvature)
@variable(model,U[1:N-1])
@variable(model,X[1:N,1:6], start = 10.0)

N = 100
s = 1:length(track.curvature)
t_final = X[end,6]


for k = 1:N-1 # loop over control intervals
    x_next = X[k,:] + (s[k+1]-s[k]) * f(s[k], U[k,:], X[k,:]);
    @constraint(model,X[k+1,:] .== x_next)
end

#initial states
#        [vx vy psi dpsi n t]
initial = [4; 0; 0; 0; 0; 0]

#set initial states
#@constraint(model,X[1,:] .== initial)

@constraint(model,X[1,1] .== 4)
@constraint(model,X[1,6] .>= 0)
@constraint(model,diff(X[:,6]) .>=0)


@objective(model,Min,X[end,6])
optimize!(model)

# Plotting the results
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Step", ylabel = "Velocity (vx)")
lines!(ax, 1:N, value.(X[:,1]))
display(fig)









