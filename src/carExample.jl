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
    
    #print(Fy)
    FxMaxsquared = max((Fz*μ + m*9.81*μ)^2 - Fy^2,0)
    # Add optimization constraints if model is provided
    if optiModel !== nothing
        @constraint(optiModel, Fy^2 + inputForce^2 <= (Fz*μ + m*9.81*μ)^2)
        FxPowerMax = maxPower/(vx)

        @constraint(optiModel,inputForce<=FxPowerMax)
        #@constraint(optiModel, inputForce <= 1000)
        #@constraint(optiModel, inputForce >= -1000)
        #@constraint(optiModel, inputForce^2 <= FxMaxsquared) #Fx = min(inputForce, FxMax)
       # @constraint(optiModel, Fx <= inputForce) #Fx = min(inputForce, FxMax)

       # Create return structure, these first 5 states are forced on all models
       #[dvx dvy dpsi ddpsi dn]
        dstates =  [inputForce/m 0  0  0 0]
    else
        Fx = min(inputForce, sqrt(FxMaxsquared))
        Fx = max(Fx, -sqrt(FxMaxsquared))
        Fx = min(Fx,maxMotorForce)
        #power limitation
        FxPowerMax = maxPower/(vx+1)
        Fx = min(Fx,FxPowerMax)
        Fx -= 1/2 *rho*CD*vx^2
#
    #    # Create return structure, these first 5 states are forced on all models
        dstates =  [Fx/m 0  0    0     0 0]
    end          # [dvx dvy dpsi ddpsi dn]
    
    

    

    return dstates

end

