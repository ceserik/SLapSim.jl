
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

inputs = car.controlMapping(car.carParameters,[1,2])

struct carStates
       speed
end

struct track
    curvature
    rho
    μ
end

trackParameters = track(
    4,
    1.225,
    1
)

states = carStates(
    23
) 

car.carFunction(car,inputs,states,trackParameters)