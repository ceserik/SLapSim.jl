include("carParams.jl")
using Interpolations



function massPointCar(car,track,k, optiModel=nothing)
    # Simple  Mass point friction circle
    # inputs is a struct of all posible inputs of vehicle, controls are not separated from inputs
    
    # Get inputs
    m          = car.carParameters.mass.value
    inputForce = car.carParameters.motorForce.value
    CL         = car.carParameters.CL.value
    CD         = car.carParameters.CD.value
    maxPower   = car.carParameters.powerLimit.value
    # Get state
    vx         = car.carParameters.vx.value
    car.carParameters.psi.value = track.theta[k]
    # Get track inputs
    c          = track.curvature[k]
    rho        = track.rho[k]
    μ          = track.μ[k]
    # Calculate forces
    Fz = 1/2 * rho * CL * vx^2
    Fy = m*vx^2*c

    maxMotorForce = 6507 #calculated for ctu25 should be added to inputs, vx torque characerisitic
    FxPowerMax = maxPower/(vx)
    #print(Fy)
    FxMaxsquared = max((Fz*μ + m*9.81*μ)^2 - Fy^2,0)
    # Add optimization constraints if model is provided
    if optiModel !== nothing
        ## tu urobit funkciu do ktorej dam obmedzenie a ona mi o spravi aby som nemusel stale davat ify
        ## aj ked to mozno je jedno lebo tento mass point bude mozno malo pouzivany?
        ## ale radsej to spravit, nech netreba prepisovat model lebo sa z toho zblaznim ked tam bude nieco inak
        @constraint(optiModel, (Fy/maxMotorForce)^2 + (inputForce/maxMotorForce)^2 <= ((Fz*μ + m*9.81*μ)/maxMotorForce)^2)
        @constraint(optiModel,inputForce/FxPowerMax <= FxPowerMax/FxPowerMax)
    else
        inputForce = min(inputForce, sqrt(FxMaxsquared))
        inputForce = max(inputForce, -sqrt(FxMaxsquared))
        inputForce = min(inputForce,maxMotorForce)
        inputForce = min(inputForce,FxPowerMax)
    end         
        dstates =  [(inputForce -(1/2 *rho*CD*vx^2))/m 0  0    0     0 0]
    return dstates
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

    function mapping(car,u,x)
        car.carParameters.motorForce.value = u[1]
        car.carParameters.vx.value = x[1]
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


    
