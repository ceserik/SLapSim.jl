using Revise
using SLapSim
using GLMakie
using CairoMakie
using Interpolations

setup_plot_theme!()
GLMakie.closeall()

# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
ipopt_attrs = Dict{String,Any}(
    "linear_solver"           => "mumps",
    "print_timing_statistics" => "yes",
    "alpha_for_y"             => "safer-min-dual-infeas",
)

results_dir = "sync/thesisImages"
mkpath(results_dir)

function run_case(max_iter::Int, suffix::String)
    path = "tracks/FSCZ.kml"
    track = kml2track(path, false, true)
    car = createTwintrack(true, track)

    exp = Experiment(
        car   = car,
        track = track,
        discipline = Open(v_start=5.0),
        solver = IpoptBackend(
            performSensitivity = false,
            attributes = ipopt_attrs,
        ),
        mesh_refinement = MeshRefinementConfig(
            tol            = 1e-1,
            method         = :h,
            error_method   = :ode,
            max_iterations = max_iter,
            segments       = Int64(round(track.sampleDistances[end] / 10)),
            pol_order      = 2,
            variant        = "Radau",
        ),
        global_constraints = GlobalConstraint[],
        analysis = AnalysisConfig(
            plot_path           = false,
            plot_controls       = false,
            plot_states         = false,
            plot_jacobian       = false,
            plot_hessian        = false,
            animate             = false,
            time_simulation     = false,
            plot_initialization = false,
            animation_path      = joinpath(results_dir, "mesh_fscz_$(suffix).mp4"),
        ),
    )

    run_experiment!(exp)

    # Sampling density
    sampling_density = get_sampling_density(exp.optiResult.segment_edges)
    fig_density, _, _ = plot_on_path(exp, sampling_density, "Mesh sampling density")
    CairoMakie.save(joinpath(results_dir, "mesh_refinement_fscz_density_$(suffix).pdf"), fig_density)
    CairoMakie.save(joinpath(results_dir, "mesh_refinement_fscz_density_$(suffix).png"), fig_density)

    # ODE residual error per segment
    errs = getSegmentErrors(exp; method=:ode)
    edges = exp.optiResult.segment_edges
    s_mid = [(edges[i] + edges[i+1]) / 2 for i in 1:length(edges)-1]
    s_knots = vcat(edges[1], s_mid, edges[end])
    err_knots = vcat(errs[1], errs, errs[end])
    node_error_itp = extrapolate(interpolate((s_knots,), err_knots, Gridded(Linear())), Flat())
    fig_error, _, _ = plot_on_path(exp, node_error_itp, "ODE residual error")
    CairoMakie.save(joinpath(results_dir, "mesh_refinement_fscz_error_$(suffix).pdf"), fig_error)
    CairoMakie.save(joinpath(results_dir, "mesh_refinement_fscz_error_$(suffix).png"), fig_error)

    return exp
end

exp_iter0  = run_case(0,  "iter0")
exp_iter10 = run_case(10, "iter10")

println("Saved mesh refinement plots to $results_dir")

nothing
