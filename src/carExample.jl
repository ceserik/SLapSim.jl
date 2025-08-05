using JuMP

function massPointCar(car, parameters, trackParameters, optiModel=nothing)
    # Simple  Mass point friction circle
    # parameters is a struct of all posible parameters of vehicle, controls are not separated from parameters
    
    # Get parameters
    m = parameters.mass.value
    inputForce = parameters.motorForce.value
    CL = parameters.CL.value
    CD = parameters.CD.value
    maxPower = parameters.powerLimit.value
    # Get state
    vx = parameters.speed.value

    # Get track parameters
    c = trackParameters.curvature
    rho = trackParameters.rho
    μ = trackParameters.μ

    # Calculate forces
    Fz = 1/2 * rho * CL * vx^2
    Fy = vx^2 * c


    maxMotorForce = 6507 #calculated for ctu25
    
    

    #print(Fy)
    FxMax = sqrt((Fz*μ + m*9.81*μ)^2 - Fy^2)
    # Add optimization constraints if model is provided
    if optiModel !== nothing
        @constraint(optiModel, inputForce <= 1000)
        @constraint(optiModel, inputForce >= -1000)
        @constraint(optiModel, Fx <= FxMax) #Fx = min(inputForce, FxMax)
        @constraint(optiModel, Fx <= inputForce) #Fx = min(inputForce, FxMax)
    
    else
        Fx = min(inputForce, FxMax)
        Fx = max(Fx, -FxMax)

        Fx = min(Fx,maxMotorForce)

        #power limitation
        FxPowerMax = maxPower/vx

        Fx = min(Fx,FxPowerMax)
        
        
        Fx -= 1/2 *rho*CD*vx^2
    end


    
    

    # Create return structure
    dstates = (ax = Fx/m,)

    return dstates
end


