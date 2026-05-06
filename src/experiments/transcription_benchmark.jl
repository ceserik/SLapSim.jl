using SLapSim
import SLapSim: solve_success
using Statistics

const TRANSCRIPTION_VARIANTS = ["Radau", "Legendre", "LobattoIIIA_integral"]

function _default_ipopt_attrs()
    return Dict{String,Any}(
        "linear_solver" => "mumps",
        "print_timing_statistics" => "yes",
        "alpha_for_y" => "safer-min-dual-infeas",
    )
end

function _default_segments(track)
    return max(1, Int(round(track.sampleDistances[end] / 10)))
end

function _build_experiment(track, car; variant::String, tol::Real, method::Symbol,
    max_iterations::Int, pol_order::Int, segments::Int, ipopt_attrs::Dict{String,Any},
    output_dir::String)

    animation_path = joinpath(output_dir, "bench_$(variant).mp4")
    analysis = AnalysisConfig(
        plot_path = false,
        plot_controls = false,
        plot_states = false,
        plot_jacobian = false,
        plot_hessian = false,
        animate = false,
        animation_speedup = 1.0,
        time_simulation = false,
        plot_initialization = false,
        animation_path = animation_path,
    )

    return Experiment(
        car = car,
        track = track,
        discipline = Open(v_start=5.0),
        solver = IpoptBackend(
            performSensitivity = false,
            attributes = ipopt_attrs,
        ),
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

function _run_case(track_name::String, track_fn, car_name::String, car_fn, variant::String;
    tol::Real, method::Symbol, max_iterations::Int, pol_order::Int, ipopt_attrs::Dict{String,Any},
    output_dir::String)

    try
        track = track_fn()
        car = car_fn(track)
        segments = _default_segments(track)
        exp = _build_experiment(track, car; variant=variant, tol=tol, method=method,
            max_iterations=max_iterations, pol_order=pol_order, segments=segments,
            ipopt_attrs=ipopt_attrs, output_dir=output_dir)

        solve_time = @elapsed run_experiment!(exp)
        status = termination_status(exp.model)
        success = solve_success(exp.model)
        obj = success ? objective_value(exp.model) : NaN

        return (
            track = track_name,
            car = car_name,
            variant = variant,
            success = success,
            status = string(status),
            solve_time = solve_time,
            objective = obj,
        )
    catch err
        return (
            track = track_name,
            car = car_name,
            variant = variant,
            success = false,
            status = "error: $(err)",
            solve_time = NaN,
            objective = NaN,
        )
    end
end

function run_transcription_benchmarks(; variants=TRANSCRIPTION_VARIANTS, tol=1e-1,
    method::Symbol=:h, max_iterations::Int=10, pol_order::Int=2,
    ipopt_attrs=_default_ipopt_attrs(), output_dir::String="results/benchmark")

    cases = [
        (
            track = "FSCZ",
            track_fn = () -> kml2track("tracks/FSCZ.kml", false, true),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
            ],
        ),
        (
            track = "Berlin",
            track_fn = () -> csv2track("src/Track/berlin_2018.csv"),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
                ("formula", t -> formulaE2026(t)),
            ],
        ),
        (
            track = "FSG",
            track_fn = () -> kml2track("tracks/FSG.kml", false, true),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
            ],
        ),
        (
            track = "DoubleTurn",
            track_fn = () -> doubleTurn(false, 0.1),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
            ],
        ),
    ]

    results = []
    for case in cases
        for (car_name, car_fn) in case.cars
            for variant in variants
                push!(results, _run_case(case.track, case.track_fn, car_name, car_fn, variant;
                    tol=tol, method=method, max_iterations=max_iterations,
                    pol_order=pol_order, ipopt_attrs=ipopt_attrs, output_dir=output_dir))
            end
        end
    end

    return results
end

function compare_transcription_methods(results)
    variants = sort(unique(getfield.(results, :variant)))
    summaries = []
    for variant in variants
        subset = [r for r in results if r.variant == variant]
        total = length(subset)
        successes = count(r -> r.success, subset)
        failures = total - successes
        times = [r.solve_time for r in subset if r.success && isfinite(r.solve_time)]

        summary = (
            variant = variant,
            runs = total,
            successes = successes,
            failures = failures,
            success_rate = total == 0 ? 0.0 : successes / total,
            mean_time = isempty(times) ? NaN : mean(times),
            median_time = isempty(times) ? NaN : median(times),
            min_time = isempty(times) ? NaN : minimum(times),
            max_time = isempty(times) ? NaN : maximum(times),
        )
        push!(summaries, summary)
    end

    return summaries
end
