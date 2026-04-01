module SLapSim
using Revise
using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
using OrdinaryDiffEq
using Debugger
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

include("animation/drawCar.jl")
include("dataAnalysis/validation.jl")
include("solvingMethods/myCollocation.jl")
include("solvingMethods/collocation.jl")


# Export public functions
export Car, Track, Result, createCTU25_1D, singleTurn, findOptimalTrajectory,kml2track, massPointSolver, createSimplestSingleTrack, time2path,initializeSolution, createTwintrack, createBasicWheelAssembly,figureEight, carParameter, CarParameters
export JuMP, Ipopt,plotCarPath, doubleTurn,plotTrack,timeSimulation, carVar, interp1, createLobattoIIIA, create_gauss_pseudospectral_metod,find_optimal_trajectory2,get_diff_matix, skidpad, create_gauss_legendre
export Drivetrain, Chassis, Motor, Gearbox, Tire, Aero, Suspension, WheelAssembly, Accumulator
export createCTU25gearbox, createFischerMotor, createR20lin, createPepikCTU25, createBasicAero, createSimpleSuspension, createDummySuspension, createCTU25chassis
export createBus, createBusSuspension, createBusWheelAssembly, RHO_SEA_LEVEL
export drawCar!, animateCar, animateCarDual, draw!
#println("SLapSim module loaded")

end # module