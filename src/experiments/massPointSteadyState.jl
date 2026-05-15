
using SLapSim
using CairoMakie

setup_plot_theme!()
# solver solves first forward pass, then bakcward pass and takes minimums of speeds
car = createCTU25_1D()
track = 0
track = singleTurn(50.0, 5.0)
track = doubleTurn()
#path = "tracks/FSCZ.kml"
#track = kml2track(path,true,false)

#smooth_by_OCP(track,0.01,0.5)
N = length(track.curvature)
track.rho = fill(track.rho[1],N)
track.μ   = fill(track.μ[1],N)


velocityfig, actualSpeed, vxMax, vForward, vBackward = massPointSolver(car,track)

results_dir = joinpath(@__DIR__, "..", "..", "sync", "massPointSteadyState")
mkpath(results_dir)

CairoMakie.save(joinpath(results_dir, "speed_profile.svg"), velocityfig)
CairoMakie.save(joinpath(results_dir, "speed_profile.pdf"), velocityfig)

trackfig = Figure(size=(500, 300))
trackax = Axis(trackfig[1, 1], aspect=DataAspect())
plotTrack(track; b_plotStartEnd=true, ax=trackax)
CairoMakie.save(joinpath(results_dir, "track.svg"), trackfig)
CairoMakie.save(joinpath(results_dir, "track.pdf"), trackfig)