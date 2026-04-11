using Revise
using SLapSim
using GLMakie
import MathOptInterface as MOI
using UnoSolver
using UnicodePlots
using DiffOpt

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

problem = Problem_config(nothing, nothing, nothing, nothing,nothing)


#car = createSimplestSingleTrack()

#car = createBus()

#track = figureEight(true, 0.1)
#track = singleTurn(50.0,5.0,true)
track = doubleTurn(true,0.1)

path = "tracks/FSCZ.kml"
#track = kml2track(path, false, true)
#track = doubleTurn(false,0.1)
#track = skidpad(false)
problem.track = track
car = createTwintrack(true,track)
#car = formulaE2026()
#car = createBus()
problem.car = car

plotTrackStates(track)
#UnoSolver.Optimizer

model = make_ipopt_model()

problem.model = model
#model = JuMP.Model(() -> UnoSolver.Optimizer(preset="ipopt"))
#optiResult = findOptimalTrajectory(track,car,model,sampleDistances,initialization)
segments = Int64(round(track.sampleDistances[end]/2))
pol_order = 2
#optiResult, optiResult_interp = find_optimal_trajectory2(problem,segments,pol_order,"Radau")
t_solve = @elapsed begin
    optiResult, optiResult_interp = find_optimal_trajectory_adaptive(problem, segments, pol_order, "Lobatto")
end
println("Solve time: $(round(t_solve, digits=2)) s")
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
#    SLapSim.plotCarStates_interp(optiResult_interp, 0.1)
    #SLapSim.plotCarStates2(optiResult)

    include("../dataAnalysis/jacobian.jl")
    include("../dataAnalysis/hessian_test2.jl")
    println("jacobian")
    println(UnicodePlots.spy(jacobian))
    println("hessian")
    println(UnicodePlots.spy(H_star))
    fig_jac = Figure()
    ax_jac = Axis(fig_jac[1, 1], title="Jacobian")
    spy!(ax_jac, SparseArrays.sparse(rotr90(jacobian)))
    display(GLMakie.Screen(), fig_jac)

    #error_itp = getErrors(problem)
    plot_controls_on_path(problem, optiResult_interp)
    sampling_density = get_sampling_density(optiResult_interp.path)
    
    plot_on_path(problem,sampling_density,"sampling density")

    snapshots = snapshot_car(car, optiResult_interp, track)
    #plot_parameters(snapshots, car,    ["drivetrain.motors[1].torque" , "drivetrain.motors[2].torque" ,"drivetrain.motors[3].torque","drivetrain.motors[4].torque"],"wheelAssemblies[1].steeringAngle")
    
    sensitivityAnalysis(problem)
    animateCarDual(track, optiResult_interp, car; speedup=1, view_radius= car.chassis.wheelbase.value*3,cam_offset=3.0, savepath="results/animation.mp4")
    plot_states_controls(car, optiResult_interp, track; sample_step=0.1)

end