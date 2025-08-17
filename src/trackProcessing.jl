using Interpolations
using JuMP, Ipopt, Zygote
using DSP
include("trackDefinition.jl")

using Interpolations


#interpolation of theta and curvature
#smart sampling by initilization from massPointSimulation

function naiveProcessing(track)
    nOfSamples = length(track.x)
    

    dx = diff(track.x)
    dy = diff(track.y)

    ds = sqrt.(dx .^ 2 .+ dy .^ 2)
    
    s = [ 0.0 ;cumsum(ds)]

    theta = atan.(dy, dx) ./ ds
    curvature = diff(theta)
    #curvature = [curvature[1]; curvature] #the curvature is shifted left by differentiating theta
    # so I copied the first element to shift it to right
    # for theta  icopied last vlaues, probably not correct but OK for now

    theta = [theta ;theta[end]]
    curvature = [curvature ; curvature[end]]

    track.theta = theta
    track.curvature = curvature
    track.sampleDistances = ds
    track.s = s
    
    return track
end

function interpolateTrack(instantTrack,track)
    th = [track.theta; track.theta[end]]
    curv = [track.curvature; track.curvature[end]]

    s = [ 0.0 ;cumsum(track.sampleDistances)]
    theta_itp = linear_interpolation(s,th)
    curv_itp  = linear_interpolation(s, curv)
end

function smooth_by_OCP(track,r,ds)

#function [x_traj, y_traj, s_traj, th_traj, C_traj] = smoothTrackByOCP(x_smpl, y_smpl, r, closedTrack, ds)
    nargin = 4
    closedTrack = false
    if nargin < 3
        r = 1e6;
    end
    
    if nargin < 4
        closedTrack = false;
    end
    
    x_smpl = track.x;
    y_smpl = track.y;
    
    dx = diff(x_smpl);
    dy = diff(y_smpl);
    ds_tmp = sqrt.(dx.^2 + dy.^2);

    s_tmp = [0; cumsum(ds_tmp)];
    
    ####
    N  = Int(round((s_tmp[end]-s_tmp[1])/ds));
    ds = (s_tmp[end]-s_tmp[1])/N;
    s_traj = collect(LinRange(s_tmp[1], s_tmp[end], Int(N)))
    x_smpl = interp1(s_tmp, x_smpl, s_traj);
    y_smpl = interp1(s_tmp, y_smpl, s_traj);    
    ###

    dx = diff(x_smpl);
    dy = diff(y_smpl);
    
    th_init = unwrap(atan.(dy, dx));
    th_init = [th_init; th_init[end]];
    C_init = diff(th_init)./ds;
    C_init =[C_init; C_init[end]];

    ####
    # Initialize the optimization problem
    model = Model(Ipopt.Optimizer) 
    
    # define decision variables



    @variable(model, u[i=1:N-1, 1], start = (diff(C_init))[i])
    @variable(model,Z[1:N,1:4])
    z_C    = Z[:, 1];
    z_th   = Z[:, 2];
    z_x    = Z[:, 3];
    z_y    = Z[:, 4];

    #u = opti.variable(N-1, 1);
    
    # Objective function (minimize the final time)
    x_dev = ((z_x[1:end-1] - x_smpl[1:end-1]).^2 + (z_x[2:end] - x_smpl[2:end]).^2)/2;
    y_dev = ((z_y[1:end-1] - y_smpl[1:end-1]).^2 + (z_y[2:end] - y_smpl[2:end]).^2)/2;


    @objective(model,Min,sum(ds .* ( r.*u.^2 .+ x_dev .+ y_dev)))
    #opti.minimize( sum(ds * ( r*u.^2 + x_dev + y_dev)) );
    
    # Dynamic constraints

    function f(z,u)
        return [u;z[1];cos(z[2]);sin(z[2])]';
    end
    #f = @(z, u) [u;z(1);cos(z(2));sin(z(2))]';

    # ---- RK4 --------
    for k=1:N-1 # loop over control intervals
        dsk = (s_traj[k+1]-s_traj[k]);
        # Runge-Kutta 4 integration
        k1 = f(Z[k, :],            u[k, :]);
        k2 = f(Z[k, :] + vec(dsk/2*k1), u[k, :]);
        k3 = f(Z[k, :] + vec(dsk/2*k2), u[k, :]);
        k4 = f(Z[k, :] + vec(dsk*k3),   u[k, :]);
        x_next = Z[k, :] + vec(dsk/6*(k1+2*k2+2*k3+k4));
        #opti.subject_to(Z[k+1, :]==x_next);# % close the gaps
        @constraint(model,Z[k+1, :]==x_next)
    end

    if closedTrack
        @constraint(model,z_x[end]  == z_x[1], start = x_smpl)
        @constraint(model,z_y[end]  == z_y[1],start = y_smpl)
        @constraint(model,z_t[end] ==z_th[1] + 2*pi*round((th_init[end]-th_init[1])/2/pi),start = th_init)
        @constraint(model,z_C[end]  == z_C[1],start = C_init)
    end
    
    #opti.set_initial(z_x, x_smpl);
    #opti.set_initial(z_y, y_smpl);
    #opti.set_initial(z_C, C_init);
    #opti.set_initial(z_th, th_init);
    
#    opti.set_initial(u, diff(C_init))
    
#     opti.subject_to(-0.05 < u < 0.05)
#     Rmax = 5;
#     opti.subject_to(-1/Rmax < z_C < 1/Rmax)
#     
#     dev_max = 2;
#     opti.subject_to( x_dev < dev_max^2 )
#     opti.subject_to( y_dev < dev_max^2 )
    


    # Solve NLP
    #opti.solver('ipopt'); % use IPOPT solver
    optimize!(model)
    #%%
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
    lines(track.x,track.y, axis = (aspect = DataAspect(),))
    return (x_traj, y_traj, C_traj, th_traj)

end



function curvateby3pointsCircle()


end


#track = singleTurn()
#smooth_by_OCP(track,1e1,0.5)

#