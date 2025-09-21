using SLapSim


car = createSimplestSingleTrack()
track = 0
track = singleTurn()
#path = "tracks/FSCZ.kml"
#track = kml2track(path,true)
model = JuMP.Model(Ipopt.Optimizer)


findOptimalTrajectory(track,car,model)