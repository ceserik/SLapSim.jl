#using Revise
using SLapSim
using Interpolations
using JuMP, Ipopt, Zygote
using DSP
using KML:KMLFile,read
using Statistics
using Proj
using PCHIPInterpolation
using Dierckx
#using LobattoInterpolation

using Interpolations

function interp1(X, V, Xq,type=nothing)
    if type == "PCHIP"
        spl = Spline1D(X, V)
        return spl(Xq)
    else
        knots = (X,)
        itp = interpolate(knots, V, Gridded(Linear()))
        itp.(Xq)
    end
end



#interpolation of theta and curvature
#smart sampling by initilization from massPointSimulation

function naiveProcessing(track::Track)
    nOfSamples = length(track.x)

    dx = diff(track.x)
    dy = diff(track.y)

    ds = sqrt.(dx .^ 2 .+ dy .^ 2)

    s = [0.0; cumsum(ds)]

    theta = atan.(dy, dx) ./ ds
    curvature = diff(theta)
    #curvature = [curvature[1]; curvature] #the curvature is shifted left by differentiating theta
    # so I copied the first element to shift it to right
    # for theta  icopied last vlaues, probably not correct but OK for now

    theta = [theta; theta[end]]
    curvature = [curvature; curvature[end]]

    track.theta = theta
    track.curvature = curvature
    track.sampleDistances = ds
    track.s = s

    return track
end

function smooth_by_OCP(track::Track, r::Float64, ds::Float64,closedTrack::Bool)
    x_smpl = track.x
    y_smpl = track.y

    dx = diff(x_smpl)
    dy = diff(y_smpl)
    ds_tmp = sqrt.(dx.^2 + dy.^2)

    s_tmp = [0; cumsum(ds_tmp)]

    ####
    N = Int(round((s_tmp[end] - s_tmp[1]) / ds))
    ds = (s_tmp[end] - s_tmp[1]) / N
    s_traj = collect(LinRange(s_tmp[1], s_tmp[end], Int(N)))
    x_smpl = interp1(s_tmp, x_smpl, s_traj,"PCHIP") 
    y_smpl = interp1(s_tmp, y_smpl, s_traj,"PCHIP") 
    ###

    dx = diff(x_smpl)
    dy = diff(y_smpl)

    th_init = unwrap(atan.(dy, dx))
    th_init = [th_init; th_init[end]]
    C_init = diff(th_init) ./ ds
    C_init = [C_init; C_init[end]]

    ####
    # Initialize the optimization problem
    model = JuMP.Model(Ipopt.Optimizer)
        function f(z::Vector{JuMP.VariableRef}, u::Vector{JuMP.VariableRef},s::Float64)
            return [u; z[1]; cos(z[2]); sin(z[2])]
        end

        function f(z::Vector{Float64}, u::Vector{Float64},s::Float64)
            return [u; z[1]; cos(z[2]); sin(z[2])]
        end

    Lobattostage = 3
    lobotom = createLobattoIIIA(Lobattostage,f)
    
    samplingPoints = length(s_traj)

    spl = Spline1D(s_traj, s_traj)
    

    xd = lobotom.createConstraints(f,4,1,spl,s_traj,model,[C_init  th_init x_smpl y_smpl],diff(C_init))
    Xall = xd[2]
    Uall = xd[3]
    Z = xd[4]
    u = xd[5]
    s_all = xd[6]
    
    z_C = Z[:, 1]
    z_th= Z[:, 2]
    z_x = Z[:, 3]
    z_y = Z[:, 4]

    # Objective function (minimize the final time)
    x_dev = ((z_x[1:end-1] - x_smpl[1:end-1]).^2 + (z_x[2:end] - x_smpl[2:end]).^2) / 2
    y_dev = ((z_y[1:end-1] - y_smpl[1:end-1]).^2 + (z_y[2:end] - y_smpl[2:end]).^2) / 2

    @objective(model, Min, sum(ds .* (r .* u.^2 .+ x_dev .+ y_dev)))

    # Dynamic constraints
    

    if closedTrack
        @constraint(model, z_x[end] == z_x[1])
        @constraint(model, z_y[end] == z_y[1])
        @constraint(model, z_th[end] == z_th[1] + 2 * pi * round((th_init[end] - th_init[1]) / 2 / pi))
        @constraint(model, z_C[end] == z_C[1])
    end

    # Solve NLP
    optimize!(model)
    # Extract solution values


    x_traj = value.(z_x)
    y_traj = value.(z_y)
    C_traj = value.(z_C)
    th_traj = value.(z_th)

    itp = lobotom.createInterpolator(value(Xall),value(Uall),s_all)
    #ssss = LinRange(s_all[1],s_all[end],546)
    
    # plot all four state components on the same axis
    #fig_interp = Figure()
    #ax_interp = Axis(fig_interp[1,1], xlabel = "s", ylabel = "State value", title = "State Interpolation - LobattoIIIA")
    #colors = (:red, :blue, :green, :black)
    #labels = ("curvature", "theta", "x", "y")
    #for stav in 1:4
    #    scatter!(ax_interp, s_all, value.(Xall[:,stav]), label = "Actual - "*labels[stav], color = colors[stav], markersize = 5)
    #    lines!(ax_interp, ssss, itp(collect(s_all[1]:0.1:s_all[end]))[:,stav], label = "Interpolated - "*labels[stav], color = colors[stav])
    #end
    #axislegend(ax_interp)
    #display(GLMakie.Screen(), fig_interp)

    track.x = x_traj
    track.y = y_traj
    track.curvature = C_traj
    track.theta = th_traj
    track.sampleDistances = s_traj    
    track.fcurve = itp
    return (x_traj, y_traj, C_traj, th_traj)

end


function plotTrackStates(track::Track)
    s = track.sampleDistances
    vals = track.fcurve.(s)              # vector of 4‑tuples: (x,y,theta,curvature)
    x  = collect(getindex.(vals, 3))
    y  = collect(getindex.(vals, 4))
    th = collect(getindex.(vals, 2))
    C  = collect(getindex.(vals, 1))

    fig = Figure(resolution = (800, 900))
    ax1 = Axis(fig[1, 1], xlabel = "s", ylabel = "x", title = "x(s)")
    ax2 = Axis(fig[2, 1], xlabel = "s", ylabel = "y", title = "y(s)")
    ax3 = Axis(fig[3, 1], xlabel = "s", ylabel = "theta", title = "theta(s)")
    ax4 = Axis(fig[4, 1], xlabel = "s", ylabel = "curvature", title = "curvature(s)")

    lines!(ax1, s, x)
    lines!(ax2, s, y)
    lines!(ax3, s, th)
    lines!(ax4, s, C)

    # show figure
    display(GLMakie.Screen(), fig)
    return fig, (ax1, ax2, ax3, ax4)
end

function kml2cart(path::String)
    file = read(path, KMLFile)
    coordinates = file.children[1].Features[1].Geometry.coordinates
    coords = hcat(ntuple(i -> getindex.(coordinates, i), 3)...)
    medLatitude = median(coords[:, 2])
    medLon = median(coords[:, 1])
    trans = Proj.Transformation("EPSG:4326", "+proj=merc +lat_ts=$(medLatitude) +lon_0=$(medLon)")
    A = zeros(Base.size(coords, 1), 3)
    for i = 1:(Base.size(coords, 1))
        points = trans(coords[i, 1], coords[i, 2])
        A[i, 1] = points[1]
        A[i, 2] = points[2]

    end
    center = trans(medLon, medLatitude)
    A[:, 1] .-= center[1]
    A[:, 2] .-= center[2]

    return A
end

function kml2track(path::String,closeTrack::Bool,flip    )
    A = kml2cart(path)

    if flip == true
        A = A[:,[2, 1, 3]]
    end


    if closeTrack
        A = [A[1:end,:]; A[1,:]']
    end
    track = Track(
        [0.0],
        [1.225],
        [1.0],
        [1.0],
        trackMapping,
        A[:,1],
        A[:,2],
        [0.0],
        [1.5],
        [1.5],
        [0.0],
        [0.0],
        (s) -> 0.0,
        [0.0]
        )
        trackfig = Figure()
        ax = Axis(trackfig[1,1],aspect=DataAspect())
        #lines!(ax,track.y, track.x)
        #display(fig)

        smooth_by_OCP(track,1.0,1.0,closeTrack)
        #track.fcurve = make_fcurve(track.sampleDistances, track.x, track.y, track.theta, track.curvature)
        #lines!(ax,track.y, track.x)
        #display(GLMakie.Screen(),trackfig)
        plotTrack(track)
    return track


end


function plotTrack(track::Track; b_plotStartEnd::Bool = false, ax::Union{Axis,Nothing} = nothing)
    # create Figure/Axis if none provided
    s = track.sampleDistances
    created = false
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1,1], aspect = DataAspect())
        created = true
    end

    # get centerline and heading from track

    vals = track.fcurve.(s)         # vector of 4-tuples
    xc  = getindex.(vals, 3)
    yc  = getindex.(vals, 4)
    thc = getindex.(vals, 2)

    # ignore 4th component or collect it similarly with getindex.(vals,4)   
    xc = collect(xc)
    yc = collect(yc)
    thc = collect(thc)

    # lane limits (broadcast to match shapes)

    xc_lim1 = - track.widthL .* sin.(thc)
    xc_lim2 = - track.widthR .* sin.(thc)

    yc_lim1 =  track.widthL .* cos.(thc)
    yc_lim2 =  track.widthR .* cos.(thc)
    ## centerline (dashed gray) and boundaries (black)
    lines!(ax, xc, yc; linestyle = :dash, linewidth = 1)
    lines!(ax, xc .+ xc_lim1, yc .+ yc_lim1; color = :black, linewidth = 1)
    lines!(ax, xc .- xc_lim2, yc .- yc_lim2; color = :black, linewidth = 1)

    ## optional start / end markers
    if b_plotStartEnd
        scatter!(ax, [xc[1]], [yc[1]]; color = :red, marker = :circle, markersize = 8)
        scatter!(ax, [xc[end]], [yc[end]]; color = :red, marker = :x, markersize = 10)
    end

    if created
        display(GLMakie.Screen(), fig)
        return fig, ax
    else
        return ax
    end
end

function make_fcurve(s_traj::Vector{Float64}, x_traj::Vector{Float64}, y_traj::Vector{Float64}, th_traj::Vector{Float64}, C_traj::Vector{Float64})
    return s -> (
        #interp1(s_traj, x_traj, s,"PCHIP"), # x position of a point
        #interp1(s_traj, y_traj, s,"PCHIP"), # y position of a point
        #interp1(s_traj, th_traj, s,"PCHIP"),# track heading at point
        #interp1(s_traj, C_traj,  s,"PCHIP") # curvate of track at a point
#
        interp1(s_traj, x_traj, s,), # x position of a point
        interp1(s_traj, y_traj, s,), # y position of a point
        interp1(s_traj, th_traj, s,),# track heading at point
        interp1(s_traj, C_traj,  s,) # curvate of track at a point
    )
end

