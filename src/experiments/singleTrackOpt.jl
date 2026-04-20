using Revise
using SLapSim
using GLMakie
using HSL_jll                      # exposes HSL_jll.libhsl_path for Ipopt's hsllib attribute
using CUDA
using MadNLP
using MadNLPGPU


#dark theme detector for linux KDE with kde-cli-tools installed
setup_plot_theme!()

GLMakie.closeall()

# ---------------------------------------------------------------------------
# Scenario setup
# ---------------------------------------------------------------------------
#car_fn = createSimplestSingleTrack
#car_fn = createBus

#track = figureEight(true, 0.1)
#track = singleTurn(50.0,5.0,true)
#track = doubleTurn(true,0.1)

path = "tracks/FSCZ.kml"
track = kml2track(path, false, true)
#track = csv2track("src/Track/berlin_2018.csv")
#track = doubleTurn(false, 0.1)
#track = skidpad(false)
#car = createSimplestSingleTrack(track)
#car = formulaE2026(track)
car = createTwintrack(true, track)
#car = formulaE2026()
#car = createBus()

plotTrackStates(track)

# ---------------------------------------------------------------------------
# Experiment configuration — flip flags to change behavior.
# ---------------------------------------------------------------------------
# See https://coin-or.github.io/Ipopt/OPTIONS.html for the full list.
ipopt_attrs = Dict{String,Any}(
    "linear_solver"          => "mumps",

    # "hsllib"                 => HSL_jll.libhsl_path,  # needed when linear_solver is ma27/57/97
    #"max_iter"               => 3000,
    #"tol"                    => 1e-2,
    #"constr_viol_tol" => 1e-2,
    #"compl_inf_tol"   => 1e-2,
    #"dual_inf_tol"    => 1e-1,
    #"mu_strategy"            => "adaptive",
    #"mu_init"                => 1e-2,
    #"acceptable_tol"         => 1e-4,
    #"acceptable_iter"        => 10,
    "print_timing_statistics"=> "yes",
    #"hessian_approximation" => "limited-memory",
    "alpha_for_y" => "safer-min-dual-infeas"
)

madnlp_attrs = Dict{String,Any}(
    "array_type"    => CUDA.CuArray,
    "linear_solver" => MadNLPGPU.CUDSSSolver,
)

exp = Experiment(
    car = car,
    track = track,
    discipline = Open(v_start=5.0),             # or Closed() for periodic BCs
    solver = IpoptBackend(
        performSensitivity = true,
        attributes = ipopt_attrs,
    ),
    #solver = MadNLPBackend(
    #    attributes = madnlp_attrs,
    #),
    # Example: cap total drive energy at 10 MJ. Leave vector empty for no global constraints.
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
        animation_path  = "sync/power.mp4",
    ),
)

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------
segments = Int64(round(track.sampleDistances[end] / 2))
pol_order = 2
run_experiment!(exp, segments, pol_order; variant="Lobatto")


sampling_density = get_sampling_density(exp.optiResult.path); # this may silently change emanigwhen higher than 2 pol order is used may cause bugs
plot_on_path(exp, sampling_density, "sampling density");

snapshots = snapshot_car(car, exp.optiResult, track)

plot_parameters(snapshots, car,
    ["drivetrain.motors[1].power", "drivetrain.motors[2].power",
     "drivetrain.motors[3].power", "drivetrain.motors[4].power"],
    ["drivetrain.accumulators.power"]
)

nothing  # suppress REPL echo of the large CarSnapshot vector
