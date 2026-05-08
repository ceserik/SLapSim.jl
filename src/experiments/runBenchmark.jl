include("transcription_benchmark.jl")

#solver_results = run_transcription_benchmarks(
#    linear_solvers=LINEAR_SOLVERS,
#    tracks=["DoubleTurn"],          # or ["FigureEight"] — small synthetic
#    cars=["singletrack"],           # cheapest car
#    max_iterations=3,               # fewer mesh refines
#)



# transcription methods (default mumps solver)
#results = run_transcription_benchmarks()
#summary = compare_transcription_methods(results)
#export_results(results, summary; output_dir="sync/thesisTables")

# linear solvers (fixed Radau variant)
solver_results = run_transcription_benchmarks(linear_solvers=LINEAR_SOLVERS)
solver_summary = compare_transcription_methods(solver_results)
export_results(solver_results, solver_summary; output_dir="sync/thesisTables/linearSolvers")






