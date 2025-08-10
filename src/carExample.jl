using JuMP

function massPointCar(car, inputs, trackParameters, optiModel=nothing)
    # Simple  Mass point friction circle
    # inputs is a struct of all posible inputs of vehicle, controls are not separated from inputs
    
    # Get inputs
    m = inputs.mass.value
    inputForce = inputs.motorForce.value
    CL = inputs.CL.value
    CD = inputs.CD.value
    maxPower = inputs.powerLimit.value
    # Get state
    vx = inputs.vx.value
    # Get track inputs
    c = trackParameters.curvature
    rho = trackParameters.rho
    μ = trackParameters.μ

    # Calculate forces
    Fz = 1/2 * rho * CL * vx^2
    Fy = vx^2 * c*m

    maxMotorForce = 6507 #calculated for ctu25 should be added to inputs, vx torque characerisitic
    FxPowerMax = maxPower/(vx)
    #print(Fy)
    FxMaxsquared = max((Fz*μ + m*9.81*μ)^2 - Fy^2,0)
    # Add optimization constraints if model is provided
    if optiModel !== nothing
        ## tu urobit funkciu do ktorej dam obmedzenie a ona mi o spravi aby som nemusel stale davat ify
        ## aj ked to mozno je jedno lebo tento mass point bude mozno malo pouzivany?
        @constraint(optiModel, Fy^2 + inputForce^2 <= (Fz*μ + m*9.81*μ)^2)
        @constraint(optiModel,inputForce<=FxPowerMax)
    else
        inputForce = min(inputForce, sqrt(FxMaxsquared))
        inputForce = max(inputForce, -sqrt(FxMaxsquared))
        inputForce = min(inputForce,maxMotorForce)
        inputForce = min(inputForce,FxPowerMax)
        
    end          # [dvx dvy dpsi ddpsi dn]

        #inputForce -= 1/2 *rho*CD*vx^2 #apply drag
        dstates =  [(inputForce -(1/2 *rho*CD*vx^2))/m 0  0    0     0 0]
    return dstates

end

