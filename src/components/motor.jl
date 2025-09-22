mutable struct Motor
    torque::carParameter{carVar}
    angularFrequency::carParameter{carVar}
    mass::carParameter{carVar}
    loss::carParameter{carVar}
    torqueSpeedFunction::Function #mapping speed to max torque
    constraints::Function
end


function createFischerMotor() 
    torque = carParameter{carVar}(0.0,"motor torque","Nm")
    angularFrequency = carParameter{carVar}(0.0,"angular frequency","rad/s")
    loss = carParameter{carVar}(0.0,"loss","W")
    torqueSpeedFunction = angularFrequency::Float64 -> 29.0
    mass = carParameter{carVar}(3.0,"mass??","kg")

    function constraints(u,optiModel::JuMP.Model)
        @constraint(optiModel, u <= 29)
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
