using SLapSim
import SLapSim: solve_success
using JuMP: termination_status, objective_value
using PrettyTables
using Printf
using Statistics

const TRANSCRIPTION_VARIANTS = ["Radau", "Legendre", "LobattoIIIA_integral"]
const LINEAR_SOLVERS = ["mumps", "ma57", "ma27", "ma97"]

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
    output_dir::String, result_label::String=variant)

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
            variant = result_label,
            success = success,
            status = string(status),
            solve_time = solve_time,
            objective = obj,
        )
    catch err
        return (
            track = track_name,
            car = car_name,
            variant = result_label,
            success = false,
            status = "error: $(err)",
            solve_time = NaN,
            objective = NaN,
        )
    end
end

function run_transcription_benchmarks(; variants=TRANSCRIPTION_VARIANTS,
    linear_solvers=nothing, tol=1e-1, method::Symbol=:h, max_iterations::Int=10,
    pol_order::Int=2, ipopt_attrs=_default_ipopt_attrs(),
    output_dir::String="results/benchmark", tracks=nothing, cars=nothing)

    cases = [
        (
            track = "FSCZ",
            track_fn = () -> kml2track("tracks/FSCZ.kml", false, true),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
#                ("bus", t -> createBus(t)),
            ],
        ),
        (
            track = "Berlin",
            track_fn = () -> csv2track("src/Track/berlin_2018.csv"; vis=false),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
                ("formula", t -> formulaE2026(t)),
                ("bus", t -> createBus(t)),
            ],
        ),
        (
            track = "FSG",
            track_fn = () -> kml2track("tracks/FSG.kml", false, true),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
#                ("bus", t -> createBus(t)),
            ],
        ),
        (
            track = "DoubleTurn",
            track_fn = () -> doubleTurn(false, 0.1),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
                ("bus", t -> createBus(t)),
            ],
        ),
        (
            track = "FigureEight",
            track_fn = () -> figureEight(false, 0.1),
            cars = [
                ("singletrack", t -> createSimplestSingleTrack(t)),
                ("twintrack", t -> createTwintrack(true, t)),
                ("bus", t -> createBus(t)),
            ],
        ),
    ]

    tracks !== nothing && (cases = filter(c -> c.track in tracks, cases))
    solver_mode = linear_solvers !== nothing
    axis = solver_mode ? linear_solvers : variants
    results = []
    for case in cases
        case_cars = cars === nothing ? case.cars : filter(c -> c[1] in cars, case.cars)
        for (car_name, car_fn) in case_cars
            for label in axis
                attrs = copy(ipopt_attrs)
                variant = solver_mode ? "Radau" : label
                solver_mode && (attrs["linear_solver"] = label)
                push!(results, _run_case(case.track, case.track_fn, car_name, car_fn, variant;
                    tol=tol, method=method, max_iterations=max_iterations,
                    pol_order=pol_order, ipopt_attrs=attrs, output_dir=output_dir,
                    result_label=label))
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

function _fmt_float(v, _i, _j)
    if v isa AbstractFloat
        isnan(v) && return "NaN"
        return @sprintf("%.4f", v)
    end
    return v
end

_short_variant(variant::AbstractString) = variant == "LobattoIIIA_integral" ? "Lobatto IIIA" : variant

function _short_status(status)
    status_str = String(status)
    startswith(status_str, "error:") && return "error"
    status_str == "LOCALLY_SOLVED" && return "solved"
    return replace(lowercase(status_str), '_' => ' ')
end

function _results_table(results)
    return (
        track = [r.track for r in results],
        car = [r.car for r in results],
        var = [_short_variant(r.variant) for r in results],
        ok = [r.success ? "yes" : "no" for r in results],
        status = [_short_status(r.status) for r in results],
        time_s = [r.solve_time for r in results],
        obj = [r.objective for r in results],
    )
end

function _summary_table(summary)
    return (
        var = [_short_variant(r.variant) for r in summary],
        runs = [r.runs for r in summary],
        ok = [r.successes for r in summary],
        fail = [r.failures for r in summary],
        rate = [100 * r.success_rate for r in summary],
        mean_s = [r.mean_time for r in summary],
        med_s = [r.median_time for r in summary],
        min_s = [r.min_time for r in summary],
        max_s = [r.max_time for r in summary],
    )
end

function _write_resized_tex(path::String, table)
    buf = IOBuffer()
    pretty_table(buf, table; backend=:latex, formatters=[_fmt_float])
    body = String(take!(buf))
    body = replace(body,
        "\\begin{tabular}" => "\\begin{longtable}",
        "\\end{tabular}" => "\\end{longtable}")
    open(path, "w") do io
        write(io, "{\\footnotesize\n")
        write(io, body)
        write(io, "}\n")
    end
end

function export_results(results, summary; output_dir::String="results/benchmark",
    show_terminal::Bool=true)

    mkpath(output_dir)
    results_table = _results_table(results)
    summary_table = _summary_table(summary)

    if show_terminal
        println("\n=== per-case results ===")
        pretty_table(results_table; formatters=[_fmt_float])
        println("\n=== variant summary ===")
        pretty_table(summary_table; formatters=[_fmt_float])
    end

    _write_resized_tex(joinpath(output_dir, "results.tex"), results_table)
    _write_resized_tex(joinpath(output_dir, "summary.tex"), summary_table)

    return joinpath(output_dir, "results.tex"), joinpath(output_dir, "summary.tex")
end
