
using SLapSim

# solver solves first forward pass, then bakcward pass and takes minimums of speeds
car = createCTU25_1D()
track = 0
#track = singleTurn()
path = "tracks/FSCZ.kml"
track = kml2track(path,true)

#smooth_by_OCP(track,0.01,0.5)
N = length(track.curvature)
track.rho = fill(track.rho[1],N)
track.μ   = fill(track.μ[1],N)


massPointSolver(car,track)