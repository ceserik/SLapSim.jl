using Revise
using SLapSim
using GLMakie
import MathOptInterface as MOI
using DiffOpt

detect = Sys.islinux()
if detect
    x = "kreadconfig6"
    option1 = "--key"
    option2 = "LookAndFeelPackage"
    try
        xd = read(`$x $option1 $option2`, String)
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

problem = Problem_config(nothing, nothing, nothing, nothing, nothing)

track = figureEight(true, 0.1)
#track = singleTurn(50.0, 5.0, true)
#track = doubleTurn(true, 0.1)
#path = "tracks/FSCZ.kml"
#track = kml2track(path, false, true)

problem.track = track
car = formulaE2026(track)
problem.car = car

plotTrackStates(track)

model = make_ipopt_model()
problem.model = model

segments = Int64(round(track.sampleDistances[end] / 2))
pol_order = 2

t_solve = @elapsed begin
    optiResult, optiResult_interp = find_optimal_trajectory_adaptive(problem, segments, pol_order, "Lobatto")
end
println("Solve time: $(round(t_solve, digits=2)) s")
problem.optiResult = optiResult_interp

if 1 == 1
    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())
    plotCarPath_interpolated(track, optiResult_interp, ax)
    screen = display(GLMakie.Screen(), fig)

    getError([1, 5], problem)
    sol = timeSimulation_interpolated(car, optiResult_interp, track)
    lines!(ax, getindex.(sol.u, 5), getindex.(sol.u, 6), label="Simulated in time")
    axislegend(ax, position=:rt)

    fig = nothing
    ax = nothing

    snapshots = snapshot_car(car, optiResult_interp, track)

    sensitivityAnalysis(problem)

    plot_parameters(snapshots, car,
        "drivetrain.motors[1].torque",
        "wheelAssemblies[1].steeringAngle",
        ["drivetrain.tires[1].brakingForce", "drivetrain.tires[2].brakingForce",
         "drivetrain.tires[3].brakingForce", "drivetrain.tires[4].brakingForce"]
    )
    plot_parameters(snapshots, car,
        "carParameters.velocity" => 1,
         "carParameters.velocity" =>2,
        "carParameters.psi",
        "carParameters.angularVelocity" => 3,
        "carParameters.n"
    )
    plot_parameters(snapshots, car,
        ["drivetrain.tires[1].forces" => 3, "drivetrain.tires[2].forces" => 3],
        ["drivetrain.tires[3].forces" => 3, "drivetrain.tires[4].forces" => 3]
    )

    animateCarDual(track, optiResult_interp, car; speedup=1,
        view_radius=car.chassis.wheelbase.value * 3, cam_offset=3.0,
        savepath="results/animation_formulaE.mp4")

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



end



