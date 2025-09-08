mutable struct WheelAssembly
    position
    velocity
    wheel2CoG
    steeringAngle::carParameter
    pivotForces::carParameter
    setPivotVelocity
    rotZ
end


function createBasicWheelAssembly(position)
    steeringAngle = carParameter(0.0,"steering angle","rad")
    pivotForces = carParameter([0.0,0.0,0.0],"Pivot forces","N N N")
    velocity = carParameter([0.0,0.0,0.0],"Velocity at pivot","")
    function rotZ()
        out = [
            cos(steeringAngle.value) -sin(steeringAngle.value) 0;
            sin(steeringAngle.value)  cos(steeringAngle.value) 0;
            0           0          1
        ]
        return out
    end

    function wheel2CoG(forces)
        #returns moments on cog from wheel forces
        moments = cross(position, forces)
        return moments
    end

    function setPivotForce(tire)
        pivotForces.value = inv(rotZ()) * tire.forces
    end

    function setPivotVelocity(angularVelocity)
        velocity.value = velocity.value + cross(angularVelocity, position) #CoG2wheelAssembly(velocity.value,angularVelocity)
    end

    function setTireSpeeds(tire)
        tire.velocity  = rotZ() * velocity
    end

    testWheelAssembly = WheelAssembly(
        position,
        [0 0 0],
        wheel2CoG,
        steeringAngle,
        pivotForces,
        setPivotVelocity,
        rotZ
        
    )
    return testWheelAssembly
end