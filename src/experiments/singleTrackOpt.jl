
using Revise
using Infiltrator
using SLapSim
using GLMakie
#Infiltrator.clear_disabled!()

#set_theme!(default_theme())

car = createSimplestSingleTrack()

#track = singleTurn(50.0,10.0,false)
path = "tracks/FSCZ.kml"
track = kml2track(path,true)


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
        #@infiltrate
        state = x[:,index]
        #
        lines!(ax,s,state,label = labels[index], linewidth = 5)
    end
    fig[1, 2] = Legend(fig, ax, "State Variables", framevisible = false)
    
    fig
end
myPlot(initialization)

model = JuMP.Model(Ipopt.Optimizer)
X,U = findOptimalTrajectory(track,car,model,initialization)
