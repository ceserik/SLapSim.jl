
using Revise
using Infiltrator
using SLapSim
using GLMakie
import MathOptInterface as MOI
using UnoSolver
using UnicodePlots
#Infiltrator.clear_disabled!()
include("../solvingMethods/optInterface.jl")
include("../dataAnalysis/validation.jl")
include("../solvingMethods/myCollocation.jl")
include("../solvingMethods/collocation.jl")
include("../solvingMethods/adaptiveRK.jl")
#include("../carModels/simpleSingleTrack.jl")
#dark theme detector for linux KDE with kde-cli-tools installed
detect = Sys.islinux()
if detect
    x = "kreadconfig6"
    option1 = "--key"
    option2 = "LookAndFeelPackage"
    try
        xd = read(`$x $option1 $option2`, String) # remember the backticks ``
        is_dark = occursin("dark", lowercase(xd))
        println("Is dark theme: ", is_dark)
        if is_dark
            set_theme!(theme_dark())
        else
            set_theme!(Makie.current_default_theme())
        end
    catch
        println("Could not detect dark theme")
    end
end

update_theme!(
    fontsize = 11,
    fonts = (; regular = "Latin Modern Roman", bold = "Latin Modern Roman Bold"),
    palette = (color = Makie.to_colormap(:tab10),),
    colormap = :turbo,
    Axis = (
        titlesize = 13,
        xlabelsize = 12,
        ylabelsize = 12,
        xticklabelsize = 10,
        yticklabelsize = 10,
    ),
    Legend = (
        labelsize = 11,
        titlesize = 12,
    ),
    Colorbar = (
        ticklabelsize = 10,
        labelsize = 11,
    ),
)

GLMakie.closeall()

problem = Problem_config(nothing, nothing, nothing, nothing)



car = createSimplestSingleTrack()
problem.car = car
track = figureEight(true, 2.0)
#track = singleTurn(50.0,5.0,true) track = doubleTurn(true,2.0)
path = "tracks/FSCZ.kml"
#htrack = kml2track(path, false, true)
#track = doubleTurn(false,0.5)
#track = skidpad(false)
problem.track = track
#Number of transcription points
#sampleDistances = collect(LinRange(track.sampleDistances[1],track.sampleDistances[end],50))
#initialization = initializeSolution(car,track,sampleDistances)
#UnoSolver.Optimizer

model = JuMP.Model(Ipopt.Optimizer)
problem.model = model
#model = JuMP.Model(() -> UnoSolver.Optimizer(preset="ipopt"))
#optiResult = findOptimalTrajectory(track,car,model,sampleDistances,initialization)
segments = 100
pol_order = 3
#optiResult, optiResult_interp = find_optimal_trajectory2(problem,segments,pol_order,"Radau")
optiResult, optiResult_interp = find_optimal_trajectory_adaptive(problem, segments, pol_order, "Radau")
problem.optiResult = optiResult_interp

if 1 == 1
    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())
    #plotCarPath(track,optiResult_interp,ax)
    plotCarPath_interpolated(track, optiResult_interp, ax)
    println(ax)
    screen = display(GLMakie.Screen(), fig)
    # simulate in time feed forward using optimal controls

    getError([1, 5], problem)
    #sol = timeSimulation(car, optiResult, track)
    sol = timeSimulation_interpolated(car, optiResult_interp, track)
    lines!(ax, getindex.(sol.u, 5), getindex.(sol.u, 6), label="Simulated in time")
    axislegend(ax, position=:rt)

    #@infiltrate
    fig = nothing
    ax = nothing
    SLapSim.plotCarStates_interp(optiResult_interp, 0.1)
    #SLapSim.plotCarStates2(optiResult)

    include("../dataAnalysis/jacobian.jl")
    include("../dataAnalysis/hessian_test2.jl")
    println("jacobian")
    println(UnicodePlots.spy(jacobian))
    println("hessian")
    println(UnicodePlots.spy(H_star))
    GLMakie.spy(rotr90(jacobian))

    error_itp = getErrors(problem)
    plot(error_itp(optiResult_interp.path))
    plotErrorsOnTrack2D(problem; itp=error_itp)

    # Plot all controls + error in a single square window (stacked rows)
    fig = Figure(resolution=(900, 900))
    # 2x2 grid of plots, each plot gets its own colorbar to the right
    plots = [
        (s -> optiResult_interp.controls(s)[1], "moment_front"),
        (s -> optiResult_interp.controls(s)[2], "moment_rear"),
        (s -> optiResult_interp.controls(s)[3], "steering"),
        (s -> error_itp(s), "error"),
    ]

    # prepare error interpolant and append as a callable



    # Layout: 2 rows × 4 columns (each plot + its colorbar occupies 2 columns)
    for i in 1:length(plots)
        row = div(i - 1, 2) + 1
        col = mod(i - 1, 2) + 1
        ax_col = 2 * col - 1
        cb_col = 2 * col
        p = plots[i]
        f = p[1]
        name = p[2]
        ax = Axis(fig[row, ax_col], aspect=DataAspect())
        _, plt = plot_on_path(problem, f, name; axis=ax)
        Colorbar(fig[row, cb_col], plt; label=name, width=25)
    end
    display(GLMakie.Screen(), fig)
    lines(error_itp(optiResult_interp.path), axis=(title="chyba v zavislosti na poloze", yscale=log10))
    sampling_density = get_sampling_density(optiResult_interp.path)
    plot_on_path(problem,sampling_density,"sampling density")


end