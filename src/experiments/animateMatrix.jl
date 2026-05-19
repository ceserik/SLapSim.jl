include("transcription_benchmark.jl")

using JuMP: termination_status
import SLapSim: solve_success

function _build_animation_experiment(track, car; tol::Real, method::Symbol,
    max_iterations::Int, pol_order::Int, segments::Int, variant::String,
    solver::SolverBackend, animation_path::String, animation_speedup::Float64)

    analysis = AnalysisConfig(
        plot_path = false,
        plot_controls = false,
        plot_states = false,
        plot_jacobian = false,
        plot_hessian = false,
        animate = true,
        animation_speedup = animation_speedup,
        time_simulation = false,
        plot_initialization = false,
        animation_path = animation_path,
    )

    return Experiment(
        car = car,
        track = track,
        discipline = Open(v_start=5.0),
        solver = solver,
        mesh_refinement = MeshRefinementConfig(
            tol = tol,
            method = method,
            error_method = :ode,
            max_iterations = max_iterations,
            segments = segments,
            pol_order = pol_order,
            variant = variant,
        ),
        global_constraints = GlobalConstraint[],
        analysis = analysis,
    )
end

"""
    run_animation_matrix(; output_dir, tol, method, max_iterations, pol_order,
                        variant, solver, animation_speedup, tracks, cars)

Solve every (track, car) combination from `_default_cases()` and write one
animation per combo into `output_dir` as `<track>_<car>.mp4`. Returns a list of
named tuples `(track, car, success, status, solve_time, animation_path)`.

Filter with `tracks=["FSCZ", ...]` and/or `cars=["singletrack", ...]`.
"""
function run_animation_matrix(; output_dir::String="sync/animations",
    tol::Real=1e-1, method::Symbol=:h, max_iterations::Int=10, pol_order::Int=2,
    variant::String="Radau", solver::SolverBackend=IpoptBackend(
        performSensitivity=false, attributes=_default_ipopt_attrs()),
    animation_speedup::Float64=1.0, tracks=nothing, cars=nothing)

    mkpath(output_dir)
    cases = _default_cases()
    tracks !== nothing && (cases = filter(c -> c.track in tracks, cases))

    results = []
    for case in cases
        case_cars = cars === nothing ? case.cars : filter(c -> c[1] in cars, case.cars)
        for (car_name, car_fn) in case_cars
            anim_path = joinpath(output_dir, "$(case.track)_$(car_name).mp4")
            entry = try
                track = case.track_fn()
                car = car_fn(track)
                segments = _default_segments(track)
                exp = _build_animation_experiment(track, car;
                    tol=tol, method=method, max_iterations=max_iterations,
                    pol_order=pol_order, segments=segments, variant=variant,
                    solver=solver, animation_path=anim_path,
                    animation_speedup=animation_speedup)

                solve_time = @elapsed run_experiment!(exp)
                status = termination_status(exp.model)
                success = solve_success(exp.model)
                (track=case.track, car=car_name, success=success,
                    status=string(status), solve_time=solve_time,
                    animation_path=anim_path)
            catch err
                (track=case.track, car=car_name, success=false,
                    status="error: $(err)", solve_time=NaN,
                    animation_path=anim_path)
            end
            push!(results, entry)
            @info "animation done" track=entry.track car=entry.car success=entry.success path=entry.animation_path
        end
    end
    return results
end



#include("src/experiments/animateMatrix.jl")

# all combos → sync/animations
results = run_animation_matrix()

# subset example:
#run_animation_matrix(tracks=["DoubleTurn", "FigureEight"], cars=["singletrack"])
