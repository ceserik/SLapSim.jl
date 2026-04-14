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
#track = doubleTurn(false, 0.1)
#track = skidpad(false)

car = createTwintrack(true, track)
#car = formulaE2026()
#car = createBus()

plotTrackStates(track)

# ---------------------------------------------------------------------------
# Experiment configuration — flip flags to change behavior.
# ---------------------------------------------------------------------------
exp = Experiment(
    car = car,
    track = track,
    discipline = Open(v_start=5.0),             # or Closed() for periodic BCs
    solver = IpoptBackend(
        linear_solver = "mumps",
        performSensitivity = false,
        print_timing = true,
        # extra_attributes = Dict("mu_init" => 1e-3),
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

# ---------------------------------------------------------------------------
# Extra analysis that isn't yet wired into run_analysis! (lives here for now).
# ---------------------------------------------------------------------------
if exp.analysis.plot_jacobian || exp.analysis.plot_hessian
    problem = exp      # legacy variable name used inside included scripts
    model = exp.model  # jacobian.jl/hessian_test2.jl reference `model` directly
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

# Sampling density plot (not yet a flag in AnalysisConfig — kept inline).
sampling_density = get_sampling_density(exp.optiResult.path);
plot_on_path(exp, sampling_density, "sampling density");

snapshots = snapshot_car(car, exp.optiResult, track)
nothing  # suppress REPL echo of the large CarSnapshot vector
