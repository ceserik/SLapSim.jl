using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
using OrdinaryDiffEq
#using Debugger 
const MODULES_INITIALIZED = true
include("carDefinition.jl")
include("massPointCar.jl")
include("simpleSingleTrack.jl")
include("massPointCar.jl")
include("trackProcessing.jl")
#include("carParams.jl")
include("aero.jl")
include("accumulator.jl")
include("gearbox.jl")
include("motor.jl")
include("tire.jl")
include("suspension.jl")
include("wheelAssembly.jl")
print("Slapsim initialized\n")


#@run car.carFunction(car)
#simplestSingleTrack(car)