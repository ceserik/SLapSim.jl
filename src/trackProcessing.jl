## Functions
# sample track from X and Y data 
# calculate vectors s, curvature, theta by OCP from the paper

#interpolation of theta and curvature


#smart sampling by initilization from massPointSimulation


function naiveProcessing(track)
    nOfSamples = length(track.x)


    dx = diff(track.x)
    dy = diff(track.y)

    ds = sqrt.(dx .^ 2 .+ dy .^ 2)

    theta = atan.(dy, dx) ./ ds
    curvature = diff(theta)
    lines(theta)

    track.theta = theta
    track.curvature = curvature


end


function smooth_by_OCP(track)

end