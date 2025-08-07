# Define the basic track structure with proper type annotations
mutable struct Track
    curvature::Any
    rho::Any
    μ::Any
    samplingDistance::Any
    mapping::Any
    x::Any
    y::Any
    theta::Any
    widthR::Any
    widthL::Any
    inclination::Any
    slope::Any
end

# Define a function to sample track properties at a given index
function trackMapping(track::Track,trackCopy::Track ,index::Integer)
    trackCopy.curvature = track.curvature[index]
end

# Constructor for a simple track with default parameters
function simpleTrack(;
    straight_length::Int = 100,
    clothoid_points::Int = 50,
    circle_length::Int = 100,
    straight_length2::Int = 100,
    ρ::Real = 1.225,  # air density
    μ::Real = 0.5,     # friction coefficient
    samplingDistance = 1.0
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
        samplingDistance,
        trackMapping,
        0,
        0,
        0
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

    angle = LinRange(pi,0,30)
    circleX = cos.(angle) .* circleRadius .+ circleRadius
    circleY = sin.(angle) .* circleRadius .+ straightLength 

    X = [straightX; circleX]
    Y = [straightY; circleY]
    #lines(X,Y)

    #track should be smoothed and curvature should be calculated properly
    #for now curvature is hardcoded

    curvature = [zeros(straightLength);fill(1/circleRadius,30)]
    theta = [fill(pi/2,straightLength);collect(angle)]
    ## will be done automatically for future tracks


    track = Track(
        curvature,
        1.225,
        1,
        1,#this is very wrong, 1 just for compatiblity, should be calculated with curvature and theta
        trackMapping,
        X,
        Y,
        theta
        )
    return track
end

# Example extension: add method to get track length
Base.length(track::Track) = length(track.curvature)
