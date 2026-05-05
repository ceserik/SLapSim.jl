using Revise
using SLapSim
using GLMakie
using HSL_jll                      # exposes HSL_jll.libhsl_path for Ipopt's hsllib attribute
using CUDA
using MadNLP
using MadNLPGPU
import MathOptInterface as MOI
using DiffOpt

setup_plot_theme!()

GLMakie.closeall()

# ---------------------------------------------------------------------------
# Scenario setup
# ---------------------------------------------------------------------------
#track = figureEight(true, 0.1)
#track = singleTurn(50.0, 5.0, true)
track = doubleTurn(true, 0.1)
#path = "tracks/FSCZ.kml"
#track = kml2track(path, false, true)
#track = csv2track("src/Track/berlin_2018.csv")

car = formulaE2026(track)

plotTrackStates(track)

# ---------------------------------------------------------------------------
# Experiment configuration — flip flags to change behavior.
# ---------------------------------------------------------------------------
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
    "alpha_for_y" => "safer-min-dual-infeas",
    #"print_level" => 8
)

madnlp_attrs = Dict{String,Any}(
    "array_type"    => CUDA.CuArray,
    "linear_solver" => MadNLPGPU.CUDSSSolver,
)

exp = Experiment(
    car = car,
    track = track,
    discipline = Open(),                        # or Closed() for periodic BCs
    solver = IpoptBackend(
        performSensitivity = true,
        attributes = ipopt_attrs,
    ),
    #solver = MadNLPBackend(
    #    attributes = madnlp_attrs,
    #),
    mesh_refinement = MeshRefinementConfig(
        tol            = 1e-1,
        method         = :h,    # :h | :p | :hp
        error_method   = :ode,  # :ode
        max_iterations = 10,
        segments       = Int64(round(track.sampleDistances[end] / 10)),
        pol_order      = 2,
        variant        = "Radau",
    ),
    global_constraints = GlobalConstraint[],    # e.g. [EnergyBudget(1.0e7)]
    analysis = AnalysisConfig(
        plot_path       = true,
        plot_controls   = true,
        plot_states     = true,
        plot_jacobian   = false,
        plot_hessian    = false,
        animate         = true,
        animation_speedup = 1.0,
        time_simulation = false,
        plot_initialization = false,
        animation_path  = "results/animation_formulaE.mp4",
    ),
)

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------
run_experiment!(exp)


sampling_density = get_sampling_density(exp.optiResult.segment_edges); # use mesh edges so density matches the sampling method
fig_sampling_density, _, _ = plot_on_path(exp, sampling_density, "sampling density")
save(joinpath(dirname(exp.analysis.animation_path), "$(splitext(basename(exp.analysis.animation_path))[1])_sampling_density.png"), fig_sampling_density)

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

nothing  # suppress REPL echo of the large CarSnapshot vector
