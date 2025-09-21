include("initSlapSim.jl")


car = createSimplestSingleTrack()
function carF(car,u,x)
    #car.mapping(car,u,x)
    car.controlMapping(car,u)
    car.stateMapping(car,x)
    dxdt = car.carFunction(car,nothing,nothing,nothing)
    return dxdt
end

#using OrdinaryDiffEq, Plots

# Example constant control function (replace with your control logic)
u_const(t) = [0.0]  # adjust length/type to match your car.controlMapping expectations

# In-place ODE wrapper expected by DifferentialEquations.jl
function f!(du, x, p, t)
    control = u_const(t)
    dx = carF(car, control, x)   # carF must return a vector matching length(x)
    #print(dx, "\n")
    du .= dx
end

# Example initial state vector: set to match what carF expects/returns
x0 = zeros(4)   # replace 6 with proper state dimension
x0[1] = 10.0    # example: initial forward speed
x0[4] = 3.0    # example: initial lateral speed
#carF(car, 4, [1 ,0, 0, 0])


tspan = (0.0, 7.0)
prob = ODEProblem(f!, x0, tspan)
sol = solve(prob, Tsit5(), saveat=0.01)

mat = hcat(sol.u...)  # size = (n_states, n_times)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "t", ylabel = "state")
for i in 1:size(mat, 1)
    lines!(ax, sol.t, vec(mat[i, :]), label = "x$(i)")
end
axislegend(ax; position = :rt)
fig