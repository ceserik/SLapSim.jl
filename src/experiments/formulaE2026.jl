using Revise
using SLapSim
using GLMakie
import MathOptInterface as MOI
using DiffOpt
using UnicodePlots
using SparseArrays

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

#track = figureEight(true, 0.1)
track = csv2track("src/Track/berlin_2018.csv")
#track = singleTurn(50.0, 5.0, true)
#track = doubleTurn(true, 0.1)
#path = "tracks/FSCZ.kml"
#track = kml2track(path, false, true)

car = formulaE2026(track)

plotTrackStates(track)

exp = Experiment(
    car = car,
    track = track,
    discipline = Closed(),                              # Formula E is a closed lap
    solver = IpoptBackend(
        linear_solver = "ma97",
        performSensitivity = false,
        print_timing = true,
    ),
    global_constraints = GlobalConstraint[],            # e.g. [EnergyBudget(1.0e7)]
    analysis = AnalysisConfig(
        plot_path       = true,
        plot_states     = true,
        plot_jacobian   = true,
        plot_hessian    = true,
        animate         = true,
        animation_path  = "results/animation_formulaE.mp4",
        time_simulation = true,
        sensitivity     = false,
    ),
)

segments = Int64(round(track.sampleDistances[end] / 6))
pol_order = 3
run_experiment!(exp, segments, pol_order; variant="Lobatto")

# ---------------------------------------------------------------------------
# Custom per-parameter plots (not yet generic enough for run_analysis!).
# ---------------------------------------------------------------------------
snapshots = snapshot_car(car, exp.optiResult, track)

plot_parameters(snapshots, car,
    "drivetrain.motors[1].torque",
    "wheelAssemblies[1].steeringAngle",
    ["drivetrain.tires[1].brakingForce", "drivetrain.tires[2].brakingForce",
     "drivetrain.tires[3].brakingForce", "drivetrain.tires[4].brakingForce"]
)
plot_parameters(snapshots, car,
    "carParameters.velocity" => 1,
    "carParameters.velocity" => 2,
    "carParameters.psi",
    "carParameters.angularVelocity" => 3,
    "carParameters.n"
)
plot_parameters(snapshots, car,
    ["drivetrain.tires[1].forces" => 3, "drivetrain.tires[2].forces" => 3],
    ["drivetrain.tires[3].forces" => 3, "drivetrain.tires[4].forces" => 3]
)

# ---------------------------------------------------------------------------
# Jacobian/Hessian inspection (inlined; included scripts reference `problem`).
# ---------------------------------------------------------------------------
if exp.analysis.plot_jacobian || exp.analysis.plot_hessian
    problem = exp                                       # legacy name used in includes
    model = exp.model                                   # jacobian.jl/hessian_test2.jl reference `model`
    include("../dataAnalysis/jacobian.jl")
    include("../dataAnalysis/hessian_test2.jl")
    if exp.analysis.plot_jacobian
        println("jacobian")
        println(UnicodePlots.spy(jacobian))
        fig_jac = Figure()
        ax_jac = Axis(fig_jac[1, 1], title="Jacobian")
        spy!(ax_jac, SparseArrays.sparse(rotr90(jacobian)))
        display(GLMakie.Screen(), fig_jac)
    end
    if exp.analysis.plot_hessian
        println("hessian")
        println(UnicodePlots.spy(H_star))
    end
end



