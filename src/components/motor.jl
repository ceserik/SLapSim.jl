mutable struct Motor{F1,F2,F3,F4}
    torque::carParameter{carVar} #actual torque
    angularFrequency::carParameter{carVar} #actual angular freuency
    mass::carParameter{carVar}
    loss::carParameter{carVar} #actual power loss
    power::carParameter{carVar} # actual mechanical power at shaft
    powerElectrical::carParameter{carVar} # electrical power drawn from battery (mech/η motoring, mech*η regen)
    efficiency::carParameter{carVar} # motor+inverter efficiency, mech<->elec
    torqueSpeedFunction::F1 #mapping speed to max torque
    constraints::F2
    setVelocity::F3
    computeElectric::F4
end
Base.show(io::IO, ::MIME"text/plain", obj::Motor) = prettyPrintComponent(io, obj)

function createFischerMotor(maxTorqueVal::Float64=29.0,maxTorqueVal_regen::Float64=29.0)
    torque = carParameter{carVar}(0.0,"motor torque","Nm",:static,[-maxTorqueVal_regen,maxTorqueVal])
    angularFrequency = carParameter{carVar}(0.0,"angular frequency","rad/s")
    loss = carParameter{carVar}(0.0,"loss","W")
    torqueSpeedFunction = angularFrequency::Float64 -> maxTorqueVal
    mass = carParameter{carVar}(3.0,"mass??","kg")
    power = carParameter{carVar}(0.0,"Power draw","W")
    powerElectrical = carParameter{carVar}(0.0,"Electrical power","W")
    efficiency = carParameter{carVar}(0.9,"motor+inverter efficiency","-",:tunable)

    function constraints(u,model=nothing)
        maxTorque = torqueSpeedFunction(0.0)
#        u = lessContraint(u/maxTorque, 1.0, model) * maxTorque
#        u = greaterContraint(u/maxTorque, -1.0, model) * maxTorque
        return u
    end

    function setVelocity(omega::carVar)
        angularFrequency.value = omega
        power.value = omega * torque.value
    end

    function computeElectric(model = nothing)
        η = efficiency.value
        p_mech = torque.value * angularFrequency.value
        if isnothing(model)
            powerElectrical.value = p_mech >= 0 ? p_mech / η : p_mech * η
        else
            # aux var pinned to max(p_mech/η, p_mech*η) by downstream battery limits
            P = @variable(model)
            @constraint(model, P >= p_mech / η)
            @constraint(model, P >= p_mech * η)
            powerElectrical.value = P
        end
    end

    motor = Motor(
        torque,
        angularFrequency,
        mass,
        loss,
        power,
        powerElectrical,
        efficiency,
        torqueSpeedFunction,
        constraints,
        setVelocity,
        computeElectric
    )
    return motor
end
