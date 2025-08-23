using Interpolations
using JuMP, Ipopt, Zygote
using DSP
using KML:KMLFile,read
using Statistics
using Proj
using PCHIPInterpolation

include("trackDefinition.jl")

using Interpolations

function interp1(X, V, Xq,type=nothing)
    

    if type == "PCHIP"
        itp = Interpolator(X, V)
        return itp.(Xq)
    else
        knots = (X,)
        itp = interpolate(knots, V, Gridded(Linear()))
        itp[Xq]
    end
end



#interpolation of theta and curvature
#smart sampling by initilization from massPointSimulation

function naiveProcessing(track)
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

function smooth_by_OCP(track, r, ds,closedTrack)
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

    # define decision variables
    @variable(model, u[i=1:N-1, 1], start = (diff(C_init))[i])
    @variable(model, Z[i =1:N,j = 1:4],start = [C_init  th_init x_smpl y_smpl][i,j])
    z_C = Z[:, 1]
    z_th= Z[:, 2]
    z_x = Z[:, 3]
    z_y = Z[:, 4]

    # Objective function (minimize the final time)
    x_dev = ((z_x[1:end-1] - x_smpl[1:end-1]).^2 + (z_x[2:end] - x_smpl[2:end]).^2) / 2
    y_dev = ((z_y[1:end-1] - y_smpl[1:end-1]).^2 + (z_y[2:end] - y_smpl[2:end]).^2) / 2

    @objective(model, Min, sum(ds .* (r .* u.^2 .+ x_dev .+ y_dev)))

    # Dynamic constraints
    function f(z, u)
        return [u; z[1]; cos(z[2]); sin(z[2])]
    end
    # ---- RK4 --------
    for k = 1:N-1 # loop over control intervals
        dsk = (s_traj[k+1] - s_traj[k])
        # Runge-Kutta 4 integration
        k1 = f(Z[k, :], u[k, :])
        k2 = f(Z[k, :] + vec(dsk / 2 * k1), u[k, :])
        k3 = f(Z[k, :] + vec(dsk / 2 * k2), u[k, :])
        k4 = f(Z[k, :] + vec(dsk * k3), u[k, :])
        x_next = Z[k, :] + vec(dsk/6*(k1 + 2*k2 + 2*k3 + k4))
        #opti.subject_to(Z[k+1, :]==x_next);# % close the gaps
        @constraint(model, Z[k+1, :] == x_next)
    end

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

    track.x = x_traj
    track.y = y_traj
    track.curvature = C_traj
    track.theta = th_traj
    track.sampleDistances = s_traj    
    return (x_traj, y_traj, C_traj, th_traj)

end

function kml2cart(path)
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

function kml2track(path,closeTrack)
    A = kml2cart(path)
    if closeTrack
        A = [A[1:end,:]; A[1,:]']
    end
    track = Track(0,1.225,1,1,trackMapping,A[:,1],A[:,2],0,0,0,0,0,0,0)
        trackfig = Figure()
        ax = Axis(trackfig[1,1],aspect=DataAspect())
        lines!(ax,track.y, track.x)
        #display(fig)

        smooth_by_OCP(track,1,1.0,closeTrack)
        lines!(ax,track.y, track.x)
        display(GLMakie.Screen(),trackfig)
    return track


end

#path = "tracks/FSCZ.kml"
#track = kml2track(path,true)
