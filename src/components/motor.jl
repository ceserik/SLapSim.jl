mutable struct Motor{F1,F2}
    torque::carParameter{carVar}
    angularFrequency::carParameter{carVar}
    mass::carParameter{carVar}
    loss::carParameter{carVar}
    torqueSpeedFunction::F1 #mapping speed to max torque
    constraints::F2
end


function createFischerMotor() 
    torque = carParameter{carVar}(0.0,"motor torque","Nm")
    angularFrequency = carParameter{carVar}(0.0,"angular frequency","rad/s")
    loss = carParameter{carVar}(0.0,"loss","W")
    torqueSpeedFunction = angularFrequency::Float64 -> 29.0
    mass = carParameter{carVar}(3.0,"mass??","kg")

    function constraints(u,model=nothing)
        maxTorque = torqueSpeedFunction(0.0)
        u = lessContraint(u/maxTorque, 29/maxTorque, model) * maxTorque
        u = greaterContraint(u/maxTorque, -29/maxTorque, model) * maxTorque
        return u
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

function createBusMotor()
    torque = carParameter{carVar}(0.0,"motor torque","Nm")
    angularFrequency = carParameter{carVar}(0.0,"angular frequency","rad/s")
    loss = carParameter{carVar}(0.0,"loss","W")
    torqueSpeedFunction = angularFrequency::Float64 -> 1000.0
    mass = carParameter{carVar}(120.0,"mass","kg")

    function constraints(u,model=nothing)
        maxTorque = torqueSpeedFunction(0.0)
        u = lessContraint(u/maxTorque, 1000/maxTorque, model) * maxTorque
        u = greaterContraint(u/maxTorque, -1000/maxTorque, model) * maxTorque
        return u
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
