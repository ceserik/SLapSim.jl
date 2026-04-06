module SLapSim
using GLMakie
using LinearAlgebra
using StaticArrays
using JuMP, Ipopt, Zygote, DiffOpt
using ControlSystemsBase
using OrdinaryDiffEq
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
include("solvingMethods/optInterface.jl")
include("dataAnalysis/sensitivityAnalysis.jl")

include("animation/drawCar.jl")
include("dataAnalysis/validation.jl")
include("solvingMethods/initialization.jl")
include("dataAnalysis/carSnapshot.jl")
include("solvingMethods/myCollocation.jl")
include("solvingMethods/collocation.jl")
include("solvingMethods/adaptiveRK.jl")

# Export public functions
export Car, Track, Result, createCTU25_1D, singleTurn, findOptimalTrajectory,kml2track, csv2track, berlinTrack, massPointSolver, createSimplestSingleTrack, time2path,initializeSolution, createTwintrack, createBasicWheelAssembly,figureEight, carParameter, CarParameters
export JuMP, Ipopt,plotCarPath, doubleTurn,plotTrack,timeSimulation, carVar, interp1, createLobattoIIIA, create_gauss_pseudospectral_metod,find_optimal_trajectory2,get_diff_matix, skidpad, create_gauss_legendre
export Drivetrain, Chassis, Motor, Gearbox, Tire, Aero, Suspension, WheelAssembly, Accumulator
export createCTU25gearbox, createFischerMotor, createR20lin, createPepikCTU25, createBasicAero, createSimpleSuspension, createDummySuspension, createCTU25chassis
export createBus, createBusSuspension, createBusWheelAssembly, RHO_SEA_LEVEL
export setParameters, resetParameters, sensitivityAnalysis, VarEntry, get_bounds, get_scales, apply_mapping!
export drawCar!, animateCar, animateCarDual, draw!
export Problem_config, Result_interpolation, find_optimal_trajectory_adaptive, make_ipopt_model,plotTrackStates
export plotCarPath_interpolated, plotCarStates_interp, getError, getErrors,get_sampling_density,plot_on_path
export timeSimulation_interpolated, snapshot_car, plot_parameters,createR20_pacejka, createQuasi_steady_Suspension,formulaE2026,plot_controls_on_path
#println("SLapSim module loaded")

end # module