using GLMakie
using LinearAlgebra
using JuMP, Ipopt, Zygote
using ControlSystemsBase
const MODULES_INITIALIZED = true
include("carCreate.jl")
include("trackProcessing.jl")
include("carParams.jl")
print("Slapsim initilized\n")

tire2= createR20lin()