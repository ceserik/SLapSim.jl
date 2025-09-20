using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
using OrdinaryDiffEq
using Debugger 
const MODULES_INITIALIZED = true
include("Track/trackProcessing.jl")
include("Track/trackDefinition.jl")
include("components/gearbox.jl")
include("components/motor.jl")
include("components/tire.jl")
include("components/accumulator.jl")
include("components/aero.jl")
include("components/suspension.jl")
include("components/wheelAssembly.jl")
include("components/chassis.jl")
include("carModels/carDefinition.jl")



include("carModels/massPointCar.jl")
include("carModels/simpleSingleTrack.jl")
include("carModels/massPointCar.jl")

#include("carParams.jl")
print("Slapsim initialized\n")

#@run car.carFunction(car)
#simplestSingleTrack(car)