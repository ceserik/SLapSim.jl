include("transcription_benchmark.jl")
using CUDA, MadNLP, MadNLPGPU

const BACKENDS = [
    "ipopt-mumps" => () -> IpoptBackend(performSensitivity=false,
        attributes=Dict{String,Any}("linear_solver"=>"mumps",
            "print_timing_statistics"=>"yes",
            "alpha_for_y"=>"safer-min-dual-infeas")),
    "madnlp-cpu" => () -> MadNLPBackend(attributes=Dict{String,Any}()),
    "madnlp-gpu" => () -> MadNLPBackend(attributes=Dict{String,Any}(
        "array_type"    => CUDA.CuArray,
        "linear_solver" => MadNLPGPU.CUDSSSolver,
        "tol"           => 1e-8,
        "acceptable_tol"=> 1e-6)),
]

#solver_results = run_transcription_benchmarks(
#    linear_solvers=LINEAR_SOLVERS,
#    tracks=["DoubleTurn"],          # or ["FigureEight"] — small synthetic
#    cars=["singletrack"],           # cheapest car
#    max_iterations=3,               # fewer mesh refines
#)



# transcription methods (default mumps solver)
results = run_transcription_benchmarks()
summary = compare_transcription_methods(results)
export_results(results, summary; output_dir="sync/thesisTables")

# linear solvers (fixed Radau variant)
solver_results = run_transcription_benchmarks(linear_solvers=LINEAR_SOLVERS)
solver_summary = compare_transcription_methods(solver_results)
export_results(solver_results, solver_summary; output_dir="sync/thesisTables/linearSolvers")

# backends: IPOPT vs MadNLP (CPU/GPU), fixed Radau variant
backend_results = run_transcription_benchmarks(backends=BACKENDS)
backend_summary = compare_transcription_methods(backend_results)
export_results(backend_results, backend_summary; output_dir="sync/thesisTables/backends")






