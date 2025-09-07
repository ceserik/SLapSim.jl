mutable struct Motor
    torque::carParameter
    angularFrequency::carParameter
    mass::carParameter
    loss::carParameter
    torqueSpeedFunction #mapping speed to max torque
    constraints
end


function createFischerMotor()
    torque = carParameter(0.0,"motor torque","Nm")
    angularFrequency = carParameter(0.0,"angular frequency","rad/s")
    loss = carParameter(0.0,"loss","W")
    torqueSpeedFunction = f(angularFrequency) = 29
    mass = carParameter(3.0,"mass??","kg")

    function constrains(u,optiModel)
        @constraint(optiModel, u <=  29)
        @constraint(optiModel, u >= -29)
    end


    motor = Motor(
        torque,
        angularFrequency,
        mass,
        loss,
        torqueSpeedFunction,
        constrains
    )
    return motor
end
