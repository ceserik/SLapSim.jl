using SLapSim
#include("../Track/trackProcessing.jl")
function skidpad(vis::Bool = false, ds::Float64 = 0.5)
        # Track parameters
        w_l = 1.5 #*0.2  # Width of the track [m]
        w_r = 1.5 # *0.2# Width of the track [m]
        rho = 1.225
        μ = 1.0

        radius = (15.25 + 1.5) / 2  # Centerline radius [m]
        straight = 10.0             # Straight section length [m]

        # Straight section before skidpad loops
        y_smpl_straight1 = collect(-straight:0)
        x_smpl_straight1 = fill(-2 * radius, length(y_smpl_straight1))

        # First loop (clockwise, two turns)
        t1 = LinRange(π, -3π, 200)
        t1 = t1[2:end]
        x_smpl_R1 = radius * cos.(t1) .- radius
        y_smpl_R1 = radius * sin.(t1)

        # Second loop mirrored to the right
        x_smpl_R2 = .-x_smpl_R1 .- 4 * radius
        y_smpl_R2 = y_smpl_R1

        # End point on top center between circles
        y_smpl_straight2 = [1.0,2.0,3.0]
        x_smpl_straight2 = [-2 * radius,-2 * radius,-2 * radius]

        # Combine all segments
        X = [x_smpl_straight1; x_smpl_R1; x_smpl_R2; x_smpl_straight2]
        Y = [y_smpl_straight1; y_smpl_R1; y_smpl_R2; y_smpl_straight2]

        track = Track(
                [0.0],                    # curvature - to be calculated
                fill(rho, length(X)),     # rho
                fill(μ, length(X)),       # μ
                [1.0],                    # sampleDistances - to be calculated properly
                trackMapping,
                X,
                Y,
                [0.0],                    # theta
                fill(w_r, length(X)),     # widthR
                fill(w_l, length(X)),     # widthL
                [0.0],                    # inclination
                [0.0],                    # slope
                s -> (0.0),               # fcurve
                [0.0]                     # s
        )

        smooth_factor = 0.5
        smooth_by_OCP(track, smooth_factor, ds, false)

        if vis == true
                plotTrack(track)
                plotTrackStates(track)
        end
        
        return track
end

#xd = skidpad(true)