
mutable struct WheelAssembly
    position::carParameter{Vector{carVar}}
    velocity::carParameter{Vector{carVar}}
    pivot2CoG::Function
    steeringAngle::carParameter{carVar}
    forces::carParameter{Vector{carVar}}
    setPivotVelocity::Function
    rotZ::Function
    setTireSpeeds::Function
    setPivotForce::Function
end


function createBasicWheelAssembly(position::Vector{carVar})
    steeringAngle = carParameter{carVar}(0.0,"steering angle","rad")
    forces = carParameter{Vector{carVar}}([0.0,0.0,0.0],"Pivot forces","N N N")
    velocity = carParameter{Vector{carVar}}([0.0,0.0,0.0],"Velocity at pivot","m/s")
    position = carParameter{Vector{carVar}}(position,"Position from CoG","m m m")
    function rotZ()
        out = [
            cos(steeringAngle.value) -sin(steeringAngle.value) 0;
            sin(steeringAngle.value)  cos(steeringAngle.value) 0;
            0           0          1
        ]
        return out
    end

    function pivot2CoG(forces)
        #returns moments on cog from wheel forces
        moments = cross(position.value, forces)
        return moments
    end

    function setPivotForce(tire::Tire)
        forces.value = inv(rotZ()) * tire.forces.value
    end

    function setPivotVelocity(angularVelocity,CoGvelocity)
        velocity.value = CoGvelocity + cross(angularVelocity, position.value) #CoG2wheelAssembly(velocity.value,angularVelocity)
    end

    function setTireSpeeds(tire::Tire)
        tire.velocity.value  = rotZ() * velocity.value
    end

    testWheelAssembly = WheelAssembly(
        position,
        velocity,
        pivot2CoG,
        steeringAngle,
        forces,
        setPivotVelocity,
        rotZ,
        setTireSpeeds,
        setPivotForce
        
    )
    return testWheelAssembly
end