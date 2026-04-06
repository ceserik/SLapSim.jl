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
    println("Smoothing Track")
    x_smpl = track.x
    y_smpl = track.y

    dx = diff(x_smpl)
    dy = diff(y_smpl)
    ds_tmp = sqrt.(dx.^2 + dy.^2)

    s_tmp = [0; cumsum(ds_tmp)]

    ####
    N = Int(round((s_tmp[end] - s_tmp[1]) / ds))
    ds = (s_tmp[end] - s_tmp[1]) / N
    s_traj = LinRange(s_tmp[1], s_tmp[end], Int(N))
    x_smpl = interp1(s_tmp, x_smpl, s_traj,"PCHIP") 
    y_smpl = interp1(s_tmp, y_smpl, s_traj,"PCHIP") 
    ###

    dx = diff(x_smpl)
    dy = diff(y_smpl)

    th_init = unwrap(atan.(dy, dx))
    th_init = [th_init[1];th_init]
    C_init = diff(th_init) ./ ds
    C_init = [C_init; C_init[end]]

    ####
    # Initialize the optimization problem
    model = JuMP.Model(Ipopt.Optimizer)
        function f(z::Vector{JuMP.VariableRef}, u::Vector{JuMP.VariableRef},s::Float64,nothing)
            return [u; z[1]; cos(z[2]); sin(z[2])]
        end

        function f(z::Union{Vector{Float64},Vector{VariableRef}}, u::Union{Vector{Float64},Vector{VariableRef}},s::Float64,nothing)
            return [u; z[1]; cos(z[2]); sin(z[2])]
        end

    Lobattostage = 2
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
    else
        @constraint(model, z_th[end] == th_init[end])
        @constraint(model, z_x[end] == x_smpl[end])
        @constraint(model, z_y[end] == y_smpl[end])
        @constraint(model, z_C[end] == C_init[end])
    end

    # Solve NLP
    #set_silent(model)
    optimize!(model)
    # Extract solution values
#    @infiltrate

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
#    @infiltrate
    return (x_traj, y_traj, C_traj, th_traj)

end


function plotTrackStates(track::Track)
    s = track.sampleDistances
    vals = track.fcurve.(s)
    C  = getindex.(vals, 1)
    th = getindex.(vals, 2)
    x  = getindex.(vals, 3)
    y  = getindex.(vals, 4)

    ds = diff(collect(s))
    u  = diff(C) ./ ds   # control: dC/ds

    fig = Figure(size = (800, 1100))
    ax1 = Axis(fig[1, 1], xlabel = "s", ylabel = "x",         title = "x(s)")
    ax2 = Axis(fig[2, 1], xlabel = "s", ylabel = "y",         title = "y(s)")
    ax3 = Axis(fig[3, 1], xlabel = "s", ylabel = "theta",     title = "theta(s)")
    ax4 = Axis(fig[4, 1], xlabel = "s", ylabel = "curvature", title = "curvature(s)")
    ax5 = Axis(fig[5, 1], xlabel = "s", ylabel = "dC/ds",     title = "control u(s)")

    lines!(ax1, s, x)
    lines!(ax2, s, y)
    lines!(ax3, s, th)
    lines!(ax4, s, C)
    lines!(ax5, collect(s)[1:end-1], u)

    display(GLMakie.Screen(), fig)
    return fig, (ax1, ax2, ax3, ax4, ax5)
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

function csv2track(path::String;
                   vis::Bool = false,
                   ds::Float64 = 2.0,
                   smooth_factor::Float64 = 1.0,
                   closedTrack::Bool = false,
                   flipXY::Bool = false,
                   widthR_default::Float64 = 1.5,
                   widthL_default::Float64 = 1.5,
                   rho::Float64 = 1.225,
                   μ::Float64 = 1.0,
                   delimiter::Char = ',',
                   comment_char::Char = '#',
                   col_x::Int = 1,
                   col_y::Int = 2,
                   col_widthR::Union{Nothing,Int} = 3,
                   col_widthL::Union{Nothing,Int} = 4)
    candidate_paths = String[path]
    if !isabspath(path)
        push!(candidate_paths, joinpath(@__DIR__, path))
    end

    csv_path = nothing
    for candidate in candidate_paths
        if isfile(candidate)
            csv_path = candidate
            break
        end
    end

    csv_path === nothing && error("Could not find track CSV. Checked: $(join(candidate_paths, ", "))")

    x = Float64[]
    y = Float64[]
    width_r = Float64[]
    width_l = Float64[]

    required_col = max(col_x, col_y)
    for line in eachline(csv_path)
        line = strip(line)
        if isempty(line) || startswith(line, string(comment_char))
            continue
        end

        cols = split(line, delimiter)
        if length(cols) < required_col
            continue
        end

        x_val = tryparse(Float64, strip(cols[col_x]))
        y_val = tryparse(Float64, strip(cols[col_y]))
        if x_val === nothing || y_val === nothing
            continue
        end

        widthR_val = widthR_default
        if col_widthR !== nothing && length(cols) >= col_widthR
            parsed_widthR = tryparse(Float64, strip(cols[col_widthR]))
            if parsed_widthR !== nothing
                widthR_val = parsed_widthR
            end
        end

        widthL_val = widthL_default
        if col_widthL !== nothing && length(cols) >= col_widthL
            parsed_widthL = tryparse(Float64, strip(cols[col_widthL]))
            if parsed_widthL !== nothing
                widthL_val = parsed_widthL
            end
        end

        push!(x, x_val)
        push!(y, y_val)
        push!(width_r, widthR_val)
        push!(width_l, widthL_val)
    end

    length(x) < 2 && error("Track CSV must contain at least two numeric (x, y) samples: $(csv_path)")

    if flipXY
        x, y = y, x
    end

    track = Track(
        [0.0],
        fill(rho, length(x)),
        fill(μ, length(x)),
        [1.0],
        trackMapping,
        x,
        y,
        [0.0],
        width_r,
        width_l,
        [0.0],
        [0.0],
        s -> 0.0,
        [0.0]
    )

    ds_raw = sqrt.(diff(x) .^ 2 .+ diff(y) .^ 2)
    s_raw = [0.0; cumsum(ds_raw)]

    smooth_by_OCP(track, smooth_factor, ds, closedTrack)

    s_new = track.sampleDistances
    track.widthR = interp1(s_raw, width_r, s_new)
    track.widthL = interp1(s_raw, width_l, s_new)
    track.rho = fill(rho, length(track.x))
    track.μ = fill(μ, length(track.x))
    track.s = s_new

    if vis
        plotTrack(track)
        plotTrackStates(track)
    end

    return track
end

function berlinTrack(; path::Union{Nothing,String} = nothing,
                     vis::Bool = false,
                     ds::Float64 = 2.0,
                     smooth_factor::Float64 = 1.0,
                     closedTrack::Bool = false,
                     flipXY::Bool = false)
    candidate_paths = String[]
    if path === nothing
        push!(candidate_paths, joinpath(@__DIR__, "berlin_208.csv"))
        push!(candidate_paths, joinpath(@__DIR__, "berlin_2018.csv"))
    else
        push!(candidate_paths, path)
        if !isabspath(path)
            push!(candidate_paths, joinpath(@__DIR__, path))
        end
    end

    csv_path = nothing
    for candidate in candidate_paths
        if isfile(candidate)
            csv_path = candidate
            break
        end
    end

    csv_path === nothing && error("Could not find Berlin track CSV. Checked: $(join(candidate_paths, ", "))")

    return csv2track(csv_path;
                     vis = vis,
                     ds = ds,
                     smooth_factor = smooth_factor,
                     closedTrack = closedTrack,
                     flipXY = flipXY)
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

        smooth_by_OCP(track,1.0,0.4,closeTrack)

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
    xc =  xc
    yc =  yc
    thc = thc

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

