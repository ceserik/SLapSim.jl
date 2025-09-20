module SLapSim

using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
using OrdinaryDiffEq
using Debugger

# Include files in dependency order
include("Track/trackProcessing.jl")
include("Track/trackDefinition.jl")

# Components (no dependencies between them)
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

# Export public functions
export Car, Track, createCTU25_1D, singleTurn, findOptimalTrajectory

println("SLapSim module loaded")

end # module