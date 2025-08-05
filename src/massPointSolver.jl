using Plots
using GLMakie
include("carCreate.jl")
include("trackDefinition.jl")
# solver solves first forward pass, then bakcward pass and takes minimums of speeds

car = createCTU25()
track = simpleTrack()


#budem sa tvarit ze disretizacia je kazdy meter
Δs = 1
vForward  = zeros(length(track.curvature))
vxMax = zeros(length(track.curvature))

#maximum speed pass
for i in eachindex(track.curvature)
    numerator = car.carParameters.mass.value * 9.81 * track.μ
    denominator = max(car.carParameters.mass.value * track.curvature[i] - 1/2 * track.rho * car.carParameters.CL.value, 0.000001)
    vxMax[i] = sqrt(numerator/denominator)  # Added sqrt for physical correctness
    vxMax[i] = min(vxMax[i], 50)  # Speed limit
end

fig  = lines(vxMax,label = "Max speed",
    axis = (; xlabel = "x", ylabel = "y", title ="Title"),
    figure = (; size = (800,600), fontsize = 22))


#create inputs for car model
inputs = deepcopy(car.carParameters)
#create instance of track parameters
trackCopy = deepcopy(track)


#forward pass
for i =  1:length(track.curvature)-1
    car.controlMapping(inputs,[99999,0])
    car.stateMapping(inputs,vForward[i])

    track.mapping(track,trackCopy,i)
    dv = car.carFunction(car,inputs,trackCopy)[1]
    dv = max(0,dv)
    v = sqrt( vForward[i]^2 +2*dv*Δs)
    vForward[i+1] = min(vxMax[i+1],v) 
    print(vForward[i+1],"\n")
    
end 



#backward pass
vBackward = vxMax
for i =  length(track.curvature):-1:2
    car.controlMapping(inputs,[-99999,0])
    car.stateMapping(inputs,vBackward[i])

    track.mapping(track,trackCopy,i)
    dv = car.carFunction(car,inputs,trackCopy)[1]
    v = sqrt(vBackward[i]^2 -2*dv*Δs)

    vBackward[i-1] = min(vxMax[i-1],v) 
    #print(v[i+1],"\n")
    
end 





lines!(vForward, label= ("vForward"))
lines!(vBackward,label = ("vBack"))
axislegend(; position = :lt)
#lines(vxMax)
#lines!(forwardSpeed)

fig


