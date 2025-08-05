# Define the basic track structure with proper type annotations
mutable struct Track{T<:Real, F<:Function}
    curvature
    rho::T
    μ::T
    mapping::F
end

# Define a function to sample track properties at a given index
function trackMapping(track::Track,trackCopy::Track ,index::Integer)
    trackCopy.curvature = track.curvature[index]
end

# Constructor for a simple track with default parameters
function simpleTrack(;
    straight_length::Int = 100,
    clothoid_points::Int = 5,
    circle_length::Int = 100,
    straight_length2::Int = 100,
    ρ::Real = 1.225,  # air density
    μ::Real = 1.0     # friction coefficient
)
    # Define track segments
    straight = fill(0.0001, straight_length)
    clothoid = LinRange(0.0, 1/15, clothoid_points)
    circle = fill(1/15, circle_length)
    straight2 = fill(0.0001, straight_length2)
    
    # Combine segments into full track
    curvature = [straight; collect(clothoid); circle;straight2]
    
    # Create and return track instance
    return Track(
        curvature,
        ρ,
        μ,
        trackMapping
    )
end

# Example extension: add method to get track length
Base.length(track::Track) = length(track.curvature)
