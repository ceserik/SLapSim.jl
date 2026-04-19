module SLapSim
using GLMakie
using CairoMakie
using LinearAlgebra
using StaticArrays
using JuMP, Ipopt, Zygote, DiffOpt
using MadNLP, MadNLPGPU
using CUDA
using HSL_jll
using ControlSystemsBase
using OrdinaryDiffEq
using DifferentiationInterface
using Interpolations
using FastGaussQuadrature
using Infiltrator
const carVar = Union{Float64, JuMP.VariableRef,JuMP.AffExpr,JuMP.NonlinearExpr,JuMP.QuadExpr}
#const carVar = Union{JuMP.VariableRef,JuMP.AffExpr,JuMP.NonlinearExpr}

# Include files in dependency order
include("Track/trackDefinition.jl")
include("Track/trackProcessing.jl")
include("Track/skidpad.jl")

# Components (no dependencies between them)
include("carModels/carParameters.jl")
include("components/gearbox.jl")
include("components/motor.jl")
include("components/tire.jl")
include("components/accumulator.jl")
include("components/aero.jl")
include("components/chassis.jl")
include("components/suspension.jl")
include("components/wheelAssembly.jl")


# Car models (depend on components)
include("carModels/carDefinition.jl")
include("carModels/massPointCar.jl")
include("carModels/simpleSingleTrack.jl")
include("carModels/multiTrack.jl")
include("carModels/bus.jl")

#solving solvingMethods
include("solvingMethods/massPointSolver.jl")
include("solvingMethods/experiment.jl")
include("solvingMethods/optInterface.jl")
include("dataAnalysis/sensitivityAnalysis.jl")

include("animation/drawCar.jl")
include("animation/themeUtils.jl")
include("dataAnalysis/validation.jl")
include("dataAnalysis/jacobian.jl")
include("dataAnalysis/hessian_test2.jl")
include("solvingMethods/initialization.jl")
include("dataAnalysis/carSnapshot.jl")
include("solvingMethods/myCollocation.jl")
include("solvingMethods/collocation.jl")
include("solvingMethods/adaptiveRK.jl")

# Export public functions
export Car, Track, Track_interpolated, interpolate_track, Result, createCTU25_1D, singleTurn, findOptimalTrajectory,kml2track, csv2track, berlinTrack, massPointSolver, createSimplestSingleTrack, time2path,initializeSolution, createTwintrack, createBasicWheelAssembly,figureEight, carParameter, CarParameters
export JuMP, Ipopt,plotCarPath, doubleTurn,plotTrack,timeSimulation, carVar, interp1, createLobattoIIIA, create_gauss_pseudospectral_metod,get_diff_matix, skidpad, create_gauss_legendre
export Drivetrain, Chassis, Motor, Gearbox, Tire, Aero, Suspension, WheelAssembly, Accumulator
export createCTU25gearbox, createFischerMotor, createR20lin, createAccumulator, createBasicAero, createSimpleSuspension, createDummySuspension, createCTU25chassis
export createBus, createBusSuspension, createBusWheelAssembly, RHO_SEA_LEVEL
export setParameters, resetParameters, sensitivityAnalysis, VarEntry, get_bounds, get_scales, apply_mapping!
export drawCar!, animateCar, animateCarDual, draw!
export detect_dark_theme!, apply_slapsim_theme!, setup_plot_theme!
export Result_interpolation, find_optimal_trajectory_adaptive, plotTrackStates
# New experiment API
export Experiment, Discipline, Open, Closed
export SolverBackend, IpoptBackend, MadNLPBackend
export AnalysisConfig, GlobalConstraint, EnergyBudget
export build_model, run_experiment!, run_analysis!
export apply_boundary_conditions!, apply_global!
export plotCarPath_interpolated, plotCarStates_interp, getError, getErrors,get_sampling_density,plot_on_path, plot_states_controls
export timeSimulation_interpolated, snapshot_car, plot_parameters,createR20_pacejka, createQuasi_steady_Suspension,formulaE2026,plot_controls_on_path
export compute_optimal_jacobian, compute_optimal_hessian, plot_jacobian_spy, plot_hessian_spy
#println("SLapSim module loaded")

end # module