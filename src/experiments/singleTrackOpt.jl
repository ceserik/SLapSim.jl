
#using Revise
using Infiltrator
using SLapSim
using GLMakie
import MathOptInterface as MOI
using UnoSolver
#Infiltrator.clear_disabled!()



#dark theme detector for linux KDE with kde-cli-tools installed
detect = true
if detect
    x = "kreadconfig6"
    option1 = "--key"
    option2 = "LookAndFeelPackage"
    xd = read(`$x $option1 $option2`,String) # remember the backticks ``
    is_dark = occursin("dark", lowercase(xd))
    println("Is dark theme: ", is_dark)
    if is_dark
        set_theme!(theme_dark())
    else
        set_theme!(Makie.current_default_theme())
    end
end


GLMakie.closeall()
car = createSimplestSingleTrack()

#track = singleTurn(50.0,5.0,true) track = doubleTurn(true,2.0)

path = "tracks/FSCZ.kml"
#track = kml2track(path,false,true)
track = doubleTurn(true,1.0)


#Number of transcription points
sampleDistances = collect(LinRange(track.sampleDistances[1],track.sampleDistances[end],10))
initialization = initializeSolution(car,track,sampleDistances)


#UnoSolver.Optimizer
model = JuMP.Model(Ipopt.Optimizer)
optiResult = findOptimalTrajectory(track,car,model,sampleDistances,initialization)


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

include("../dataAnalysis/jacobian.jl")
GLMakie.spy(reverse(jacobian,dims=1))