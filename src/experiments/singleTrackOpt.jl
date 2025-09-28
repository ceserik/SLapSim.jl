using SLapSim
using GLMakie
using JuMP

car = createSimplestSingleTrack()
track = 0
track = singleTurn(25.0,0.0)
#path = "tracks/FSCZ.kml"
#track = kml2track(path,true)
model = JuMP.Model(Ipopt.Optimizer)


X,U = findOptimalTrajectory(track,car,model)