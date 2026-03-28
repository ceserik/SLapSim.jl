using Revise
using SLapSim


include("../solvingMethods/optInterface.jl")
include("../dataAnalysis/validation.jl")
include("../solvingMethods/myCollocation.jl")
include("../solvingMethods/collocation.jl")
include("../solvingMethods/adaptiveRK.jl")
include("../carModels/multiTrack.jl")
car = createTwintrack()

#track = figureEight(true, 2.0)
car.carFunction(track)