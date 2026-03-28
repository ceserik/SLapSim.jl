using GLMakie

fig = Figure()
ax = Axis(fig[1,1], xlabel="Path", ylabel="State Value")

path_values = collect(optiResult_interp.path[1]:0.01:optiResult_interp.path[end])
states = [optiResult_interp.states(s) for s in path_values]

lines!(ax, path_values, getindex.(states, 1), label="State 1")
axislegend(ax)

display(fig)