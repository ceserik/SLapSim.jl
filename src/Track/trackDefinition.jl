using Revise
using GLMakie
using SLapSim
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

function singleTurn(straightLength::Float64,circleLength::Float64,vis::Bool = false)
    # X is East
    # Y is North
    # Z is Up
    # 0 degree around Z axis is X axis
    # ENU

    #straightLength = 5
    straightLength = Int(straightLength)
    circleRadius = 10
    straightX = zeros(straightLength)
    straightY = LinRange(0,straightLength,straightLength)

    clothoidLength = Int(circleLength)
    clothoidAngle = LinRange(pi/2, -pi/2*0, clothoidLength)
    clothoidX = zeros(clothoidLength)
    clothoidY = zeros(clothoidLength)

    for i in 1:clothoidLength
        s = i / clothoidLength
        clothoidX[i] = circleRadius * s * cos(clothoidAngle[i])
        clothoidY[i] = straightLength + circleRadius * s * sin(clothoidAngle[i])
    end


    X = [straightX; clothoidX]
    Y = [straightY; clothoidY]
    track = Track(
        [0.0],#curvature,
        [1.225],
        [1.0],
        [1.0],#this is very wrong, 1 just for compatiblity, should be calculated with curvature and theta
        trackMapping,
        X,
        Y,
        [0.0],
        [2.0],
        [2.0],
        [0.0],
        [0.0],
        s->(0.0),
        [0.0]
        )

    smooth_by_OCP(track,1.0,0.5,false)
    track.fcurve = make_fcurve(track.sampleDistances, track.x, track.y, track.theta, track.curvature)
    if vis ==1
        plotTrack(track)
    end
    return track
end

function doubleTurn(vis::Bool = false,ds::Float64 =0.5)
    # Track parameters
    w_l = 3.0  # Width of the track [m]
    w_r = 3.0  # Width of the track [m]
    
    R1 = 6.0   # radius of the first turn [m]
    R2 = 5.0   # radius of the second turn [m]
    l_straight = 20.0  # length of the straight at the beginning of the track [m]
    
    # First turn (semicircle from pi to 0)
    t1 = LinRange(π, 0, 100)
    x_smpl_R1 = R1 * cos.(t1) .- R1
    y_smpl_R1 = R1 * sin.(t1)
    
    # Second turn (semicircle from pi to 2*pi, excluding first point)
    t2 = LinRange(π, 2π, 100)
    x_smpl_R2 = R2 * cos.(t2[2:end]) .+ R2
    y_smpl_R2 = R2 * sin.(t2[2:end])
    
    # Straight section
    y_smpl_straight1 = collect(-l_straight:-1)
    x_smpl_straight1 = fill(-2*R1, length(y_smpl_straight1))
    
    # Combine all segments
    X = [x_smpl_straight1; x_smpl_R1; x_smpl_R2]
    Y = [y_smpl_straight1; y_smpl_R1; y_smpl_R2]
    
    track = Track(
        [0.0],      # curvature - to be calculated
        [1.225],    # rho
        [1.0],      # μ
        [1.0],      # sampleDistances - to be calculated properly
        trackMapping,
        X,
        Y,
        [0.0],      # theta
        fill(w_r, length(X)),  # widthR
        fill(w_l, length(X)),  # widthL
        [0.0],      # inclination
        [0.0],      # slope
        s->(0.0),   # fcurve
        [0.0]       # s
    )
    
    smooth_factor = 1e1
    smooth_by_OCP(track, 1.0,ds, false)
    
    # Update the width arrays to match the new track length after smoothing
    track.widthR = fill(w_r, length(track.x))
    track.widthL = fill(w_l, length(track.x))
    
    track.fcurve = make_fcurve(track.sampleDistances, track.x, track.y, track.theta, track.curvature)
    
    if vis == 1
        plotTrack(track)
    end
    
    return track
end

# Example extension: add method to get track length
Base.length(track::Track) = length(track.curvature)




