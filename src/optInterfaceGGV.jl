using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
include("carCreate.jl")
include("trackProcessing.jl")

function f(k,u,x)
    car.mapping(car,u,x)
    dzds = time2path(car,track,k) #time2path(s,instantCarParams,track,car)
    return dzds
end

#poradie car,track,k,s,u,x

## here I need to define trnasofmation of ODE with respect to time to ODE with respect to path
function time2path(car,track,k)
    #track.mapping(track,instantTrack,s)
    dxdt = car.carFunction(car,track,k,model) ### what the fuck preco mam instant car params a car, ked to je to iste pouzivat konzistentne poradie vo volani funkcii
    v_x = car.carParameters.vx.value
    v_y = car.carParameters.vy.value
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

car = createCTU25()
track = 0
track = singleTurn()
smooth_by_OCP(track,1e1,0.5)
N = length(track.curvature)
#fill track parameters which are constant along track, should be automatized

track.rho = fill(track.rho,N)
track.μ   = fill(track.μ,N)


#the track is assumed to be discretized before

model = Model(Ipopt.Optimizer) 
#create inputs for car model #create instance of track parameters
instantCarParams = deepcopy(car.carParameters)
instantTrack = deepcopy(track)

cas = 1000
tfinal_0 = cas
#determine sizes of inputs and states

@variable(model,U[1:N-1],start = 10.0)
@variable(model,X[1:N,1:6], start = 10.0)

N = length(track.sampleDistances)
s = track.sampleDistances#1:length(track.curvature)
t_final = X[end,6]


for k = 1:N-2 # loop over control intervals
    x_next = X[k,:] + (s[k+1]-s[k]) * f(k, U[k,:], X[k,:]);
    @constraint(model,X[k+1,:] .== x_next)
end

#initial states
#        [vx vy psi dpsi n t]
#initial = [4; 0; 0; 0; 0; 0]
#set initial states
#@constraint(model,X[1,:] .== initial)

@constraint(model,X[1,1] .== 1)
@constraint(model,X[1,6] .>= 0)
@constraint(model,diff(X[:,6]) .>=0)


@objective(model,Min,X[end,6])
optimize!(model)

# Plotting the results
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], xlabel = "Step", ylabel = "Velocity (vx)")
lines!(ax, value.(X[:,1]))
display(GLMakie.Screen(), fig)  # This creates a new window









