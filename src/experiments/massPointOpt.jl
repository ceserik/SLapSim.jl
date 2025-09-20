include("../initSlapSim.jl")
include("../solvingMethods/optInterface.jl")
car = createCTU25_1D()
track = 0
#track = singleTurn()
path = "tracks/FSCZ.kml"
track = kml2track(path,true)
model = JuMP.Model(Ipopt.Optimizer)


findOptimalTrajectory(track,car,model)