using SLapSim
using OrdinaryDiffEq
using GLMakie

car = createSimplestSingleTrack();
function carF(car,u,x)
    #car.mapping(car,u,x)
    car.controlMapping(car,u);
    car.stateMapping(car,x);
    dxdt = car.carFunction(car,nothing,nothing,nothing);
    return dxdt
end;

#using OrdinaryDiffEq, Plots

# Example constant control function (replace with your control logic)
u_const(t) = [1.0,1.0,0.2]  # adjust length/type to match your car.controlMapping expectations

# In-place ODE wrapper expected by DifferentialEquations.jl
function f!(du, x, p, t)
    car = p;
    control = u_const(t);
    dx = carF(car, control, x)  ; # carF must return a vector matching length(x)

    du .= dx;
end;


x0 = zeros(4)  ; # replace 6 with proper state dimension
x0[1] = 5  ;  # example: initial forward speed
x0[4] = 0;    # example: initial forward speed

tspan = (0.0, 10);
prob = ODEProblem(f!, x0, tspan,car);
sol = solve(prob, Tsit5(), saveat=0.01);




mat = hcat(sol.u...)  # size = (n_states, n_times)

# Define state labels
state_labels = ["vx", "vy", "psi", "dpsi"]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "t [s]", ylabel = "states")
for i in 1:size(mat, 1)
    lines!(ax, sol.t, vec(mat[i, :]), label = state_labels[i])
end
axislegend(ax; position = :rt)
fig
