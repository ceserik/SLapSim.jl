using Revise
using SLapSim
using GLMakie
using Interpolations
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

track = doubleTurn(false, 0.1)
car = createTwintrack(true, track)
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

exp = Experiment(
    car = car,
    track = track,
    discipline = Open(v_start=5.0),             # or Closed() for periodic BCs
    solver = IpoptBackend(
        performSensitivity = false,
        attributes = ipopt_attrs,
    ),

    mesh_refinement = MeshRefinementConfig(
        tol            = 1e-2,
        method         = :p,    # :h | :p | :hp
        error_method   = :ode,  # :ode
        max_iterations = 1,
        segments       = 1,#Int64(round(track.sampleDistances[end] / 10)),
        pol_order      = 20,
        variant        = "Radau",
    ),
    # Example: cap total drive energy at 10 MJ. Leave vector empty for no global constraints.
    global_constraints = GlobalConstraint[],    # e.g. [EnergyBudget(1.0e7)]
    analysis = AnalysisConfig(
        plot_path       = false,
        plot_controls   = true,
        plot_states     = true,
        plot_jacobian   = false,
        plot_hessian    = false,
        animate         = false,
        animation_speedup = 1.0,
        time_simulation = false,
        plot_initialization = false,
        animation_path  = "sync/p_methodplots/power.mp4",
    ),
)

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------
run_experiment!(exp)


#sampling_density = get_sampling_density(exp.optiResult.segment_edges); # use mesh edges so density matches the sampling method
errs, s_mid = getNodeErrors(exp; iteration=1, savepath=joinpath(dirname(exp.analysis.animation_path), "$(splitext(basename(exp.analysis.animation_path))[1])_node_errors_scatter.png"))
node_error_itp = extrapolate(interpolate((s_mid,), errs, Gridded(Linear())), Flat())
fig_node_error, _, _ = plot_on_path(exp, node_error_itp, "Node error"; colorscale=log10)
save(joinpath(dirname(exp.analysis.animation_path), "$(splitext(basename(exp.analysis.animation_path))[1])_node_error.png"), fig_node_error)
#plot_on_path(exp, sampling_density, "sampling density");

# 2x2 summary: node errors scatter | steering | torque rear | node errors on path
let
    s_path = exp.optiResult.path
    steering_itp(s)    = exp.optiResult.controls(s)[2]
    torque_rear_itp(s) = exp.optiResult.controls(s)[4]

    # Thesis sizing: ctuthesis 11pt textwidth = 363pt. Figure overflows margins
    # at 1.4*textwidth = 508pt physical in PDF. Match Makie width in pt and keep
    # fontsize = 11 so on-graph text renders at body font size after \includegraphics.
    fig_w = 508
    fig_h = round(Int, fig_w * 0.65)

    fig = Figure(size=(fig_w, fig_h), fontsize=11, figure_padding=4)
    ax_err   = Axis(fig[1, 1]; yscale=log10, title="Node-interval errors", xlabel="s [m]", ylabel="error")
    ax_path  = Axis(fig[2, 1]; aspect=DataAspect(), title="Node error on path")
    ax_steer = Axis(fig[1, 2]; aspect=DataAspect(), title="Steering angle on path")
    ax_torq  = Axis(fig[2, 2]; aspect=DataAspect(), title="Torque rear on path")

    scatter!(ax_err, s_mid, max.(errs, 1e-12))
    _, plt_err   = plot_on_path(exp, node_error_itp,  "Node error";       axis=ax_path,  colorscale=log10)
    _, plt_steer = plot_on_path(exp, steering_itp,    "Steering [rad]";   axis=ax_steer)
    _, plt_torq  = plot_on_path(exp, torque_rear_itp, "Brake [F]"; axis=ax_torq)

    Colorbar(fig[2, 0], plt_err;   label="Node error",       width=8, flipaxis=false)
    Colorbar(fig[1, 3], plt_steer; label="Steering [rad]",   width=8)
    Colorbar(fig[2, 3], plt_torq;  label="Brake [F]", width=8)

    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 4)

    display(GLMakie.Screen(), fig)
    save(joinpath(dirname(exp.analysis.animation_path), "$(splitext(basename(exp.analysis.animation_path))[1])_summary2x2.png"), fig; px_per_unit=4)
end



snapshots = snapshot_car(car, exp.optiResult, track)

plot_parameters(snapshots, car,
    ["drivetrain.motors[1].power", "drivetrain.motors[2].power",
     "drivetrain.motors[3].power", "drivetrain.motors[4].power"],
    ["drivetrain.accumulators.power"]
)

nothing  # suppress REPL echo of the large CarSnapshot vector
