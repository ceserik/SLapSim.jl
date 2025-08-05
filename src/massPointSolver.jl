using Plots
using GLMakie
include("carCreate.jl")
# solver solves first forward pass, then bakcward pass and takes minimums of speeds

car = createCTU25()


# straight line then constant radius attached by clothoid
straight = zeros(100).+0.0001
clothoidCurve = LinRange(0,1/15,50)
circle = fill(1/15,100)


road = [ straight; clothoidCurve; circle]
mutable struct track
    curvature
    rho
    μ
end

trackParameters = track(
    4,
    1.225,
    1
)



#budem sa tvarit ze disretizacia je kazdy meter
Δs = 1
forwardSpeed  = zeros(length(road))
vxMax = zeros(length(road))

#maximum speed pass
for i = 1:length(road)
    numerator = car.carParameters.mass.value*9.81*trackParameters.μ
    denominator = max(car.carParameters.mass.value*road[i].-1/2*trackParameters.rho*car.carParameters.CL.value,0.000001)
    vxMax[i] = min(numerator/denominator,500)
end

fig  = lines(vxMax,label = "Max speed",
    axis = (; xlabel = "x", ylabel = "y", title ="Title"),
    figure = (; size = (800,600), fontsize = 22))

#save("./images/parabola.png", fig)



#forward pass
for i =  1:length(road)-1
    inputs = car.controlMapping(car.carParameters,[99999,0])
    inputs = car.stateMapping(inputs,forwardSpeed[i])

    trackParameters.curvature = road[i]
    dv = car.carFunction(car,inputs,trackParameters)[1]
    Δt = 2*Δs / dv
    forwardSpeed[i+1] = min(vxMax[i+1],forwardSpeed[i] + dv*Δt) 
    print(forwardSpeed[i+1],"\n")
    
end

lines!(forwardSpeed, label= ("forwardPassSpeed"))
axislegend(; position = :lt)
#lines(vxMax)
#lines!(forwardSpeed)

fig


