include("transcription_benchmark.jl")

results = run_transcription_benchmarks()
summary = compare_transcription_methods(results)
export_results(results, summary)