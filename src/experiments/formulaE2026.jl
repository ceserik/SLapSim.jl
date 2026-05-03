using Revise
using SLapSim
using GLMakie
import MathOptInterface as MOI
using DiffOpt

setup_plot_theme!()

GLMakie.closeall()

#track = figureEight(true, 0.1)
track = csv2track("src/Track/berlin_2018.csv")
#track = singleTurn(50.0, 5.0, true)
#track = doubleTurn(true, 0.1)
#path = "tracks/FSCZ.kml"
#track = kml2track(path, false, true)

car = formulaE2026(track)

plotTrackStates(track)


# See https://coin-or.github.io/Ipopt/OPTIONS.html for the full list.
ipopt_attrs = Dict{String,Any}(
    "linear_solver"          => "mumps",

    # "hsllib"                 => HSL_jll.libhsl_path,  # needed when linear_solver is ma27/57/97
    #"max_iter"               => 3000,
    #"tol"                    => 1e-6,
    #"mu_strategy"            => "adaptive",
    #"mu_init"                => 1e-2,
    #"acceptable_tol"         => 1e-4,
    #"acceptable_iter"        => 10,
    "print_timing_statistics"=> "yes",
    #"hessian_approximation" => "limited-memory",
    "alpha_for_y" => "safer-min-dual-infeas"
)

exp = Experiment(
    car = car,
    track = track,
    discipline = Open(),                              # Formula E is a closed lap
    solver = IpoptBackend(
        performSensitivity = true,
        attributes = ipopt_attrs,
    ),
    global_constraints = GlobalConstraint[],            # e.g. [EnergyBudget(1.0e7)]
    mesh_refinement = MeshRefinementConfig(
        segments  = Int64(round(track.sampleDistances[end] / 6)),
        pol_order = 3,
        variant   = "Lobatto",
    ),
    analysis = AnalysisConfig(
        plot_path       = true,
        plot_states     = true,
        plot_jacobian   = true,
        plot_hessian    = true,
        animate         = false,
        animation_path  = "results/animation_formulaE.mp4",
        time_simulation = true,
    ),
)

run_experiment!(exp)

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





