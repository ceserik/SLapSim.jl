using GLMakie
using JuMP
include("carCreate.jl")
#include("trackDefinition.jl")
include("trackProcessing.jl")
# solver solves first forward pass, then bakcward pass and takes minimums of speeds
car = createCTU25()
track = 0
track = singleTurn()
path = "tracks/FSG.kml"
track = kml2track(path,true)

#smooth_by_OCP(track,0.01,0.5)
N = length(track.curvature)
track.rho = fill(track.rho,N)
track.μ   = fill(track.μ,N)

function massPointSolver(car, track)
    Δs = diff(track.sampleDistances)
    vForward = zeros(length(track.curvature))
    vxMax = zeros(length(track.curvature))
    #create inputs for car model #create instance of track parameters
    inputs = deepcopy(car.carParameters)
    trackCopy = deepcopy(track)

    #maximum speed pass
    for i in eachindex(track.curvature)
        numerator = car.carParameters.mass.value * 9.81 * track.μ[i]
        denominator = max(car.carParameters.mass.value * abs(track.curvature[i]) - 1 / 2 * track.rho[i] * car.carParameters.CL.value, 0.000001)
        #a = track.μ[i]/car.carParameters.mass.value/abs(track.curvature[i])

        #numerator = car.carParameters.mass.value *9.81*a
        #denominator = max(1 - a*0.5*track.rho[i]*car.carParameters.CL.value,0.00001)
        vxMax[i] = sqrt(numerator / denominator) 
        vxMax[i] = min(vxMax[i], 50)  # Speed limit
    end

    #forward pass
    for i = 1:length(track.curvature)-1
        car.mapping(car,[99999],vForward[i])
        dv = car.carFunction(car,trackCopy,i)[1]
        dv = max(0, dv)
        v = sqrt(vForward[i]^2 + 2 * dv * Δs[i])
        vForward[i+1] = min(vxMax[i+1], v)
        print(vForward[i+1], "\n")

    end

    #backward pass
    vBackward = deepcopy(vxMax)
    for i = length(track.curvature)-1:-1:2
        car.mapping(car,[-99999],vBackward[i])
        dv = car.carFunction(car,trackCopy,i)[1]
        v = sqrt(vBackward[i]^2 - 2 * dv * Δs[i])

        vBackward[i-1] = min(vxMax[i-1], v)
    end


    ## plotting
    velocityfig = Figure(size=(800, 600))
    ax1 = Axis(velocityfig[1, 1],
        xlabel="Distance [m]",
        ylabel="Velocity [m/s]",
        title="Velocity Profile"
    )

    lines!(ax1, track.sampleDistances, vxMax, label="Max speed")
    lines!(ax1, track.sampleDistances, vForward, label="Forward")
    lines!(ax1, track.sampleDistances, vBackward, label="Backward")

    axislegend(ax1, position=:lt)
    display(GLMakie.Screen(), velocityfig)  # Force new window

end


massPointSolver(car,track)