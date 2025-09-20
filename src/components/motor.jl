mutable struct Motor
    torque::carParameter{Float64}
    angularFrequency::carParameter{Float64}
    mass::carParameter{Float64}
    loss::carParameter{Float64}
    torqueSpeedFunction::Function #mapping speed to max torque
    constraints::Function
end


function createFischerMotor()
    torque = carParameter(0.0,"motor torque","Nm")
    angularFrequency = carParameter(0.0,"angular frequency","rad/s")
    loss = carParameter(0.0,"loss","W")
    torqueSpeedFunction = f(angularFrequency) = 29.0
    mass = carParameter(3.0,"mass??","kg")

    function constraints(u::Union{VariableRef,Float64},optiModel)
        @constraint(optiModel, u <=  29)
        @constraint(optiModel, u >= -29)
    end


    motor = Motor(
        torque,
        angularFrequency,
        mass,
        loss,
        torqueSpeedFunction,
        constraints
    )
    return motor
end
