include("carExample.jl")
include("carParams.jl")


using Interpolations
function interp1(X, V, Xq)
    knots = (X,)
    itp = interpolate(knots, V, Gridded(Linear()))
    itp[Xq]
end

#Simple mass point model
function createCTU25()
    mass = carParameter(280.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    vx = carParameter(15.0,"Speed X","m/s")
    powerLimit = carParameter(80000.0,"PowerLimit","W")
    vy = carParameter(0.0,"Speed Y","m/s")
    psi = carParameter(0.0,"heading","rad")
    n = carParameter(0.0,"Distance from centerline","m")

    p = carParameters(
        mass,
        motorForce,
        CL,
        CD,
        vx,
        vy,
        psi,
        n,
        powerLimit
    )

    function controlMapping(input, controls)
        input.motorForce.value = controls[1]  
        #paramcopy.CL.value = controls[2]
        return input
    end

    function mapping(carParams,track,trackCopy,controls,states,s)
        carParams.motorForce.value = controls[1]
        carParams.vx.value = states[1]

        ## have to interpolate these makea function like in martinG laptimer
        #interp1(track.samplingDistance,track.theta,s)
        carParams.psi.value = interp1(track.samplingDistance,track.theta,s)
        trackCopy.curvature = interp1(track.samplingDistance,track.curvature,s)
        trackCopy.theta     = interp1(track.samplingDistance,track.theta,s)
        ## interpolate these two
    end

    function stateMapping(input,states)
        input.vx.value = states[1]
        #map psi to track heading here, maybe input should be also track mapping should be unified
        return input
    end
    car = Car(
        massPointCar,
        p,
        controlMapping,
        stateMapping,
        mapping
    )

    return car
end


    
