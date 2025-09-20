using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
using OrdinaryDiffEq
#using Debugger 
const MODULES_INITIALIZED = true
include("carModels/carDefinition.jl")
include("carModels/massPointCar.jl")
include("carModels/simpleSingleTrack.jl")
include("carModels/massPointCar.jl")
include("Track/trackProcessing.jl")
#include("carParams.jl")
include("components/aero.jl")
include("components/accumulator.jl")
include("components/gearbox.jl")
include("components/motor.jl")
include("components/tire.jl")
include("components/suspension.jl")
include("components/wheelAssembly.jl")
print("Slapsim initialized\n")


#@run car.carFunction(car)
#simplestSingleTrack(car)