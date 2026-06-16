# SLapSim.jl a Scalabe Laptime Simulator
## Modular vehicle models, ALL vehicle parameters can be chosen to be control variable, fixed parameter or optimizable parameter

One simulator for multiple complexities of vehicle models

[Thesis (PDF)](https://github.com/ceserik/SLapSim.jl/blob/main/thesis/thesis.pdf)

# Examples

[Youtube playlist](https://youtube.com/playlist?list=PL_TqVFNheNEL-B8Fof0HaNX4ybwSa8F1H)

## 2026 Class Project: Race Car Minimum Lap Time 
[Demo video (MP4)](https://github.com/user-attachments/assets/8b247832-738e-479a-9686-4933a349607b)

modelled according to https://smdogroup.github.io/ae6310/jupyter/2026_project/race_car_problem.html

## Twintrack AWD, Formula Student Germany
[Demo video (MP4)](https://github.com/user-attachments/assets/04ffdac6-0c06-49ac-b160-7d4860dca80a)


## Bus with front and rear steering
[Demo video (MP4)](https://github.com/user-attachments/assets/27b3ef73-6ede-4d06-b786-d44f40c6076a)


## Singletrack FWD
[Demo video (MP4)](https://github.com/user-attachments/assets/1c7c6e15-870a-4ea0-b19b-fa1ffae2fc7d)

This is a demo of front wheel drive vehicle on Formula Student Czech track

# Modeling framework
Power bond graphs inspired by Brown, Forbes T. Engineering System Dynamics: A Unified Graph-Centered Approach. This approach is best visualised on the singletrack model.

<img width="607" height="842" alt="image" src="https://github.com/user-attachments/assets/b28e3d06-13bd-4510-914f-efcb937c0b06" /> 


# Getting Started

### Requirements
- Julia 1.10+
- Optional: HSL solvers (`HSL_jll`) for faster Ipopt, CUDA + `MadNLPGPU` for GPU.

### Install
Clone and instantiate dependencies:
```bash
git clone https://github.com/ceserik/SLapSim.jl.git
cd SLapSim.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run an example
From the repo root:
```bash
julia --project=. src/experiments/runGenericExperiment.jl
```

Other ready-to-run scripts in [`src/experiments/`](src/experiments/):
- `massPointOpt.jl` — mass-point lap optimisation
- `runGenericExperiment.jl` — generic experiment runner (pick car model + track)
- `formulaE2026.jl` — twin-track Formula Electric on Berlin track
- `renderVehicles.jl` — animate a saved solution

### Pick a track and car
Inside any experiment script, swap the track:
```julia
track = doubleTurn(false, 0.1)       # synthetic
track = kml2track("tracks/FSCZ.kml", false, true)  # from KML
track = csv2track("src/Track/berlin_2018.csv")     # from CSV
```
and the car model:
```julia
car = createSimplestSingleTrack(track)
car = createBus(track)
```

### Output
Solver prints lap time; GLMakie window shows trajectory and states. Animations save to `sync/animations/`.

