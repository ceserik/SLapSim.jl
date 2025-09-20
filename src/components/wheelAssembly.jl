mutable struct WheelAssembly
    position::carParameter{Vector{Float64}}
    velocity::carParameter{Vector{Float64}}
    pivot2CoG::Function
    steeringAngle::carParameter{Float64}
    forces::carParameter{Vector{Float64}}
    setPivotVelocity::Function
    rotZ::Function
    setTireSpeeds::Function
    setPivotForce::Function
end


function createBasicWheelAssembly(position)
    steeringAngle = carParameter(0.0,"steering angle","rad")
    forces = carParameter([0.0,0.0,0.0],"Pivot forces","N N N")
    velocity = carParameter([0.0,0.0,0.0],"Velocity at pivot","m/s")
    position = carParameter(position,"Position from CoG","m m m")
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

    function setPivotForce(tire)
        forces.value = inv(rotZ()) * tire.forces.value
    end

    function setPivotVelocity(angularVelocity,CoGvelocity)
        velocity.value = CoGvelocity + cross(angularVelocity, position.value) #CoG2wheelAssembly(velocity.value,angularVelocity)
    end

    function setTireSpeeds(tire)
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