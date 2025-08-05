using JuMP

function massPointCar(car, parameters, states, trackParameters, optiModel=nothing)
    # Simple  Mass point friction circle
    # parameters is a struct of all posible parameters of vehicle, controls are not separated from parameters
    
    # Get parameters
    m = parameters.mass.value
    FxMotor = parameters.motorForce.value
    CL = parameters.CL.value
    
    # Get state
    vx = states.speed

    # Get track parameters
    c = trackParameters.curvature
    rho = trackParameters.rho
    μ = trackParameters.μ

    # Calculate forces
    Fz = 1/2 * rho * CL * vx^2
    Fy = vx^2 * c
    
    # Add optimization constraints if model is provided
    if optiModel !== nothing
        @constraint(optiModel, FxMotor <= 10)
        @constraint(optiModel, FxMotor >= -10)
    end
    
    FxMax = sqrt((Fz*μ + m*9.81*μ)^2 - Fy^2)
    Fx = min(FxMotor, FxMax)

    # Create return structure
    dstates = (ax = Fx/m,)

    return dstates
end


