module SLapSim

using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
using OrdinaryDiffEq
using Debugger
using Revise
using Infiltrator

const carVar = Union{Float64, JuMP.VariableRef,JuMP.AffExpr,JuMP.NonlinearExpr}
#const carVar = Union{JuMP.VariableRef,JuMP.AffExpr,JuMP.NonlinearExpr}

# Include files in dependency order
include("Track/trackDefinition.jl")
include("Track/trackProcessing.jl")

# Components (no dependencies between them)
include("carModels/carParameters.jl")
include("components/gearbox.jl")
include("components/motor.jl")
include("components/tire.jl")
include("components/accumulator.jl")
include("components/aero.jl")
include("components/suspension.jl")
include("components/wheelAssembly.jl")
include("components/chassis.jl")

# Car models (depend on components)
include("carModels/carDefinition.jl")
include("carModels/massPointCar.jl")
include("carModels/simpleSingleTrack.jl")

#solving solvingMethods
include("solvingMethods/massPointSolver.jl")
include("solvingMethods/optInterface.jl")


# Export public functions
export Car, Track, createCTU25_1D, singleTurn, findOptimalTrajectory,kml2track, massPointSolver, createSimplestSingleTrack, time2path,initializeSolution
export JuMP, Ipopt  # Re-export for convenience
#println("SLapSim module loaded")

end # module