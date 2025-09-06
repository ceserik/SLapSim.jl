mutable struct Motor
    torque::carParameter
    angularFrequency::carParameter
    mass::carParameter
    loss::carParameter
    torqueSpeedFunction #mapping speed to max torque
end


function createFischerMotor()
    torque = carParameter(0.0,"motor torque","Nm")
    angularFrequency = carParameter(0.0,"angular frequency","rad/s")
    loss = carParameter(0.0,"loss","W")
    torqueSpeedFunction = f(angularFrequency) = 29
    mass = carParameter(3.0,"mass??","kg")
    motor = Motor(
        torque,
        angularFrequency,
        mass,
        loss,
        torqueSpeedFunction
    )
    return motor
end
