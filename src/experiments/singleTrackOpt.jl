
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
initialization = initializeSolution(car,track)


model = JuMP.Model(Ipopt.Optimizer)

optiResult = findOptimalTrajectory(track,car,model,initialization)


fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
plotCarPath(track,optiResult,ax)
println(ax)
display(fig)
# simulate in time feed forward using optimal controls
sol = timeSimulation(car, optiResult, track)
lines!(ax,getindex.(sol.u, 5), getindex.(sol.u, 6))

display(fig)

fig = nothing
ax = nothing
SLapSim.plotCarStates(optiResult)
SLapSim.plotCarStates2(optiResult)