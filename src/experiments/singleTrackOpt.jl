using SLapSim


car = createSimplestSingleTrack()
track = 0
track = singleTurn(30.0)
#path = "tracks/FSCZ.kml"
#track = kml2track(path,true)
model = JuMP.Model(Ipopt.Optimizer)


findOptimalTrajectory(track,car,model)