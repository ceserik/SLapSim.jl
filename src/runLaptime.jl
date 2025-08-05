
include("carCreate.jl")
using Revise
using SLapSim
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase


car = createCTU25()
model = Model(Ipopt.Optimizer)
#car.mass.value = @variable(model,u)
#out = car.controlMapping(car.carParameters,[@variable(model,u1),@variable(model,u2)])



inputs = car.controlMapping(car.carParameters,[@variable(model,u),2])
inputs = car.stateMapping(inputs,[2])
dx = car.carFunction(car,inputs,trackParameters)