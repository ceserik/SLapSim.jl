
using Revise
using Infiltrator
using SLapSim
using GLMakie
#Infiltrator.clear_disabled!()

#set_theme!(theme_light())
GLMakie.closeall()
car = createSimplestSingleTrack()

#track = singleTurn(50.0,5.0,true)
track = doubleTurn(true)
path = "tracks/FSCZ.kml"
#track = kml2track(path,true)
#@infiltrate
result = initializeSolution(car,track)

#lines(sol)



function myPlot(result)
    x = result.states
    u = result.controls
    s = result.path
    fig = Figure()
    ax = fig[1, 1] = Axis(fig)
    labels = ["Vx", "Vy", "ψ", "ψ̇", "n", "t"] 
    #@infiltrate
    for index = 1:length(x[1,:])
        state = x[:,index]
        lines!(ax,s,state,label = labels[index], linewidth = 5)
    end
    axislegend(ax, "State Variables", position = :rt)
    fig
    fig2 = Figure()
    ax2 = fig2[1, 1] = Axis(fig2)
    ax2 = [Axis(fig2[i, 1]) for i in 1:3]
    controls = ["MomentFront", "MomentRear", "Steering"] 
    for index = 1:length(u[1,:])
        control = u[:,index]
        lines!(ax2[index], s[1:end-1], control, label = controls[index], linewidth = 5)
        axislegend(ax2[index], position = :rt)
    end
    fig2
 


end
#myPlot(result)

model = JuMP.Model(Ipopt.Optimizer)
result = findOptimalTrajectory(track,car,model,result)

myPlot(result)
plotCarPath(track,result)

myPlot(result)