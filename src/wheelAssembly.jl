mutable struct WheelAssembly
    position
    velocity
    CoG2wheelAssembly
    wheel2CoG
    steeringAngle::carParameter
    rotZ
end



function createBasicWheelAssembly(position)
    steeringAngle = carParameter(0.0,"steering angle","rad")
    function rotZ(angle)
        out = [
            cos(angle) -sin(angle) 0;
            sin(angle)  cos(angle) 0;
            0           0          1
        ]
        return out
    end

    function CoG2wheelAssembly(velocity, angularVelocity)
        #returns speed at pivot point
        return velocity + cross(angularVelocity, position)
    end
    function wheel2CoG(forces, position)
        #returns moments on cog from wheel forces
        moments = cross(position, forces)
        return moments
    end

    testWheelAssembly = WheelAssembly(
        position,
        [0 0 0],
        CoG2wheelAssembly,
        wheel2CoG,
        steeringAngle,
        rotZ
    )
    return testWheelAssembly
end