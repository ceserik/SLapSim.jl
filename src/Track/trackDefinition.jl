using GLMakie
#include("trackProcessing.jl")
# Define the basic track structure with proper type annotations
mutable struct Track
    curvature::Vector{Float64}
    rho::Vector{Float64}
    μ::Vector{Float64}
    sampleDistances::Vector{Float64} #rename to sample distances samplingDistance
    mapping::Function
    x::Vector{Float64}
    y::Vector{Float64}
    theta::Vector{Float64}
    widthR::Vector{Float64}
    widthL::Vector{Float64}
    inclination::Vector{Float64}
    slope::Vector{Float64}
    fcurve::Function
    s::Vector{Float64}
end


# interpolate track parameters at given distance s
function trackMapping(track::Track,trackCopy::Track ,s)
    trackCopy.curvature = interp1(track.sampleDistances,track.curvature,s)
    trackCopy.theta = interp1(track.sampleDistances,track.theta,s)
end

# Constructor for a simple track with default parameters
function simpleTrack(;
    straight_length::Int = 100,
    clothoid_points::Int = 50,
    circle_length::Int = 100,
    straight_length2::Int = 100,
    ρ::Real = 1.225,  # air density
    μ::Real = 0.5,     # friction coefficient
    sampleDistances = 1.0
)
    # Define track segments
    straight = fill(0.0001, straight_length)
    clothoid = LinRange(0.0, 1/10, clothoid_points)
    circle = fill(1/10, circle_length)
    straight2 = fill(0.0001, straight_length2)
    
    # Combine segments into full track
    curvature = [straight; collect(clothoid); circle;straight2]
    
    # Create and return track instance
    return Track(
        curvature,
        ρ,
        μ,
        sampleDistances,
        trackMapping,
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        s->(0.0),
        [0.0]
    )
end

function singleTurn()
    # X is East
    # Y is North
    # Z is Up
    # 0 degree around Z axis is X axis
    # ENU

    straightLength = 70
    circleRadius = 15
    straightX = zeros(straightLength)
    straightY = LinRange(0,straightLength,straightLength)

    clothoidLength = 70
    clothoidAngle = LinRange(pi/2, 0, clothoidLength)
    clothoidX = zeros(clothoidLength)
    clothoidY = zeros(clothoidLength)

    for i in 1:clothoidLength
        s = i / clothoidLength
        clothoidX[i] = circleRadius * s * cos(clothoidAngle[i])
        clothoidY[i] = straightLength + circleRadius * s * sin(clothoidAngle[i])
    end


    X = [straightX; clothoidX]
    Y = [straightY; clothoidY]

    #angle = LinRange(pi,0,30)
    #circleX = cos.(angle) .* circleRadius .+ circleRadius
    #circleY = sin.(angle) .* circleRadius .+ straightLength 

    #X = [straightX; circleX[1:end]]
    #Y = [straightY.-1; circleY[1:end]]
    #lines(X,Y)

    #track should be smoothed and curvature should be calculated properly
    #for now curvature is hardcoded

    #curvature = [zeros(straightLength);fill(1/circleRadius,30)]
    #theta = [fill(pi/2,straightLength);collect(angle)]
    ## will be done automatically for future tracks


    track = Track(
        [0.0],#curvature,
        [1.225],
        [1.0],
        [1.0],#this is very wrong, 1 just for compatiblity, should be calculated with curvature and theta
        trackMapping,
        X,
        Y,
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        s->(0.0),
        [0.0]
        )

    smooth_by_OCP(track,1,1.0,false)
    track.fcurve = make_fcurve(track.sampleDistances, track.x, track.y, track.theta, track.curvature)
    plotTrack(track,track.sampleDistances)
    return track
end

# Example extension: add method to get track length
Base.length(track::Track) = length(track.curvature)
