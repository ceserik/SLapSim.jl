using Revise
using SLapSim
using GLMakie
import MathOptInterface as MOI
using UnoSolver
using UnicodePlots
using DiffOpt
using HSL_jll                      # exposes HSL_jll.libhsl_path for Ipopt's hsllib attribute


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

# ---------------------------------------------------------------------------
# Scenario setup
# ---------------------------------------------------------------------------
#car_fn = createSimplestSingleTrack
#car_fn = createBus

#track = figureEight(true, 0.1)
#track = singleTurn(50.0,5.0,true)
#track = doubleTurn(true,0.1)

path = "tracks/FSCZ.kml"
#track = kml2track(path, false, true)
track = doubleTurn(false, 0.1)
#track = skidpad(false)

car = createTwintrack(true, track)
#car = formulaE2026()
#car = createBus()

plotTrackStates(track)

# ---------------------------------------------------------------------------
# Experiment configuration — flip flags to change behavior.
# ---------------------------------------------------------------------------
# All Ipopt options in one dict — add/remove keys freely, no struct edits needed.
# See https://coin-or.github.io/Ipopt/OPTIONS.html for the full list.
ipopt_attrs = Dict{String,Any}(
    "linear_solver"          => "ma97",
    "hsllib"                 => HSL_jll.libhsl_path,   # needed for ma27/57/97
    "max_iter"               => 3000,
    "tol"                    => 1e-6,
    "mu_strategy"            => "adaptive",
    "mu_init"                => 1e-2,
    "acceptable_tol"         => 1e-4,
    "acceptable_iter"        => 10,
    "print_timing_statistics"=> "yes",
    # "hessian_approximation" => "limited-memory",
)

exp = Experiment(
    car = car,
    track = track,
    discipline = Open(v_start=5.0),             # or Closed() for periodic BCs
    solver = IpoptBackend(
        performSensitivity = false,
        attributes = ipopt_attrs,
    ),
    # Example: cap total drive energy at 10 MJ. Leave vector empty for no global constraints.
    global_constraints = GlobalConstraint[],    # e.g. [EnergyBudget(1.0e7)]
    analysis = AnalysisConfig(
        plot_path       = true,
        plot_controls   = true,
        plot_states     = true,
        plot_jacobian   = true,
        plot_hessian    = true,
        animate         = true,
        animation_speedup = 1.0,
        time_simulation = true,
        sensitivity     = false,
    ),
)

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------
segments = Int64(round(track.sampleDistances[end] / 2))
pol_order = 2
run_experiment!(exp, segments, pol_order; variant="Lobatto")

# Jacobian/Hessian spy plots are handled by run_analysis! via
# AnalysisConfig.plot_jacobian / plot_hessian.

# Sampling density plot (not yet a flag in AnalysisConfig — kept inline).
sampling_density = get_sampling_density(exp.optiResult.path);
plot_on_path(exp, sampling_density, "sampling density");

snapshots = snapshot_car(car, exp.optiResult, track)
nothing  # suppress REPL echo of the large CarSnapshot vector
