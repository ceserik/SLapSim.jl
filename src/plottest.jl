x = 1:0.1:10
fig = lines(x, x.^2; label = "Parabola",
    axis = (; xlabel = "x", ylabel = "y", title ="Title"),
    figure = (; size = (800,600), fontsize = 22))
axislegend(; position = :lt)
#save("./images/parabola.png", fig)
fig