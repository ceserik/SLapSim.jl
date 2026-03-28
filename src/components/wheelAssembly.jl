using Revise
mutable struct WheelAssembly{F1,F2,F3}
    position::carParameter{Vector{carVar}}
    velocityPivot::carParameter{Vector{carVar}}
    velocityTire::carParameter{Vector{carVar}}
    maxAngle::carParameter{carVar}
    steeringAngle::carParameter{carVar}
    forces::carParameter{Vector{carVar}}
    torque::carParameter{Vector{carVar}}
    constraints::F1
    setVelocity::F2
    getTorque::F3
end
Base.show(io::IO, ::MIME"text/plain", obj::WheelAssembly) = prettyPrintComponent(io, obj)

function createBasicWheelAssembly(position::Vector{carVar})
    steeringAngle = carParameter{carVar}(0.0, "steering angle", "rad")
    maxAngle = carParameter{carVar}(20 / 180 * pi, " max steering angle", "rad")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Pivot forces", "N N N")
    torque = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Generated torque on CoG", "N N N")
    velocityPivot = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Velocity at pivot", "m/s")
    velocityTire = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Velocity in tire frame", "m/s")
    position = carParameter{Vector{carVar}}(position, "Position from CoG", "m m m")
    function rotZ(steering::carVar) #check direction
        out = [
            cos(steering) -sin(steering) 0;
            sin(steering) cos(steering) 0;
            0 0 1
        ]
        return out
    end

    function rotZinv(steering::carVar) #check direction
        out = [
            cos(steering) sin(steering) 0;
            -sin(steering) cos(steering) 0;
            0 0 1
        ]
        return out
    end
    function constraints(model=nothing)
        lessContraint(steeringAngle.value / maxAngle.value, 20 / 180 * pi / maxAngle.value, model)
        greaterContraint(steeringAngle.value / maxAngle.value, -20 / 180 * pi / maxAngle.value, model)
    end

    function pivot2CoG(forces)
        #returns moments on cog from wheel forces
        torque.value = cross(position.value, forces)
    end

    function setPivotForce(forcesIn)
        forces.value = rotZ(steeringAngle.value) * forcesIn
    end

    function setPivotVelocity(angularVelocity, CoGvelocity)
        velocityPivot.value = CoGvelocity + cross(angularVelocity, position.value) #CoG2wheelAssembly(velocity.value,angularVelocity)
    end

    function setTireSpeeds()
        velocityTire.value = rotZinv(steeringAngle.value) * velocityPivot.value
    end

    function setVelocity(angularVelocity, velocity)
        setPivotVelocity(angularVelocity, velocity)
        setTireSpeeds()
    end

    function getTorque(forcesIn)
        setPivotForce(forcesIn)
        pivot2CoG(forces.value)
    end

    testWheelAssembly = WheelAssembly(
        position,
        velocityPivot,
        velocityTire,
        maxAngle,
        steeringAngle,
        forces,
        torque,
        constraints,
        setVelocity,
        getTorque
    )
    return testWheelAssembly
end