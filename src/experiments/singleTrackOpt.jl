
using Revise
using Infiltrator
using SLapSim
using GLMakie
#Infiltrator.clear_disabled!()

set_theme!(theme_dark())

car = createSimplestSingleTrack()

track = singleTurn(50.0,10.0,false)
path = "tracks/FSCZ.kml"
#track = kml2track(path,true)


#@infiltrate
initialization = initializeSolution(car,track)

#lines(sol)



function myPlot(initialization)
    x = initialization.states
    u = initialization.controls
    s = initialization.path
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
    controls = ["MomentFront", "MomentRear", "Steering"] 
    for index = 1:length(u[1,:])
        control = u[:,index]
        #@infiltrate
        lines!(ax2,s[1:end-1],control,label = controls[index], linewidth = 5)
    end
    axislegend(ax2, "State Variables", position = :rt)
    fig2



end
myPlot(initialization)

model = JuMP.Model(Ipopt.Optimizer)
out = findOptimalTrajectory(track,car,model,initialization)

myPlot(out)
