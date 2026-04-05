mutable struct Motor{F1,F2}
    torque::carParameter{carVar} #actual torque
    angularFrequency::carParameter{carVar} #actual angular freuency
    mass::carParameter{carVar} 
    loss::carParameter{carVar} #actual power loss
    power::carParameter{carVar} # actual power draw
    torqueSpeedFunction::F1 #mapping speed to max torque
    constraints::F2
end
Base.show(io::IO, ::MIME"text/plain", obj::Motor) = prettyPrintComponent(io, obj)

function createFischerMotor(maxTorqueVal::Float64=29.0)
    torque = carParameter{carVar}(0.0,"motor torque","Nm",:static,[-maxTorqueVal,maxTorqueVal])
    angularFrequency = carParameter{carVar}(0.0,"angular frequency","rad/s")
    loss = carParameter{carVar}(0.0,"loss","W")
    torqueSpeedFunction = angularFrequency::Float64 -> maxTorqueVal
    mass = carParameter{carVar}(3.0,"mass??","kg")
    power = carParameter{carVar}(3.0,"Power draw","kW")

    function constraints(u,model=nothing)
        maxTorque = torqueSpeedFunction(0.0)
#        u = lessContraint(u/maxTorque, 1.0, model) * maxTorque
#        u = greaterContraint(u/maxTorque, -1.0, model) * maxTorque
        return u
    end


    motor = Motor(
        torque,
        angularFrequency,
        mass,
        loss,
        power,
        torqueSpeedFunction,
        constraints
    )
    return motor
end
