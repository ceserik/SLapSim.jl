
#using Revise
using Infiltrator
using SLapSim
using GLMakie
import MathOptInterface as MOI
using UnoSolver
using UnicodePlots
#Infiltrator.clear_disabled!()
include("../solvingMethods/optInterface.jl")
include("../dataAnalysis/validation.jl")

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

problem = Problem_config(0,0,0,0)


car = createSimplestSingleTrack()
problem.car = car
#track = singleTurn(50.0,5.0,true) track = doubleTurn(true,2.0)

path = "tracks/FSCZ.kml"
#track = kml2track(path,false,true)
track = doubleTurn(false,0.5)
problem.track = track

#Number of transcription points
#sampleDistances = collect(LinRange(track.sampleDistances[1],track.sampleDistances[end],50))
#initialization = initializeSolution(car,track,sampleDistances)



#UnoSolver.Optimizer
model = JuMP.Model(Ipopt.Optimizer)
problem.model = model
#model = JuMP.Model(() -> UnoSolver.Optimizer(preset="ipopt"))
#optiResult = findOptimalTrajectory(track,car,model,sampleDistances,initialization)
segments = 250
pol_order = 2
optiResult, optiResult_interp = find_optimal_trajectory2(track,car,model,segments,pol_order)

problem.optiResult = optiResult_interp

fig = Figure()
ax = Axis(fig[1,1], aspect = DataAspect())
#plotCarPath(track,optiResult_interp,ax)
plotCarPath_interpolated(track,optiResult_interp,ax)
println(ax)
screen = display(GLMakie.Screen(), fig)
# simulate in time feed forward using optimal controls

getError([1,5],problem)
#sol = timeSimulation(car, optiResult, track)
sol = timeSimulation_interpolated(car, optiResult_interp, track)
lines!(ax,getindex.(sol.u, 5), getindex.(sol.u, 6), label = "Simulated in time")
axislegend(ax, position = :rt)

#@infiltrate
fig = nothing
ax = nothing
SLapSim.plotCarStates_interp(optiResult_interp,0.01)
#SLapSim.plotCarStates2(optiResult)

include("../dataAnalysis/jacobian.jl")
include("../dataAnalysis/hessian_test2.jl")

println(UnicodePlots.spy(jacobian))

println(UnicodePlots.spy(H_star))
GLMakie.spy(rotr90(jacobian)) 

plot(getErrors(problem))