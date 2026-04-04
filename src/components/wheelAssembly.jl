mutable struct WheelAssembly{F1,F2,F3,F4,F5}
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
    setupObservables::F4
    updateObservables::F5
end
Base.show(io::IO, ::MIME"text/plain", obj::WheelAssembly) = prettyPrintComponent(io, obj)

function createBasicWheelAssembly(position::Vector{carVar})
    steeringAngle = carParameter{carVar}(0.0, "steering angle", "rad",:control,[-25 / 180 * pi, 25 / 180 * pi])
    maxAngle = carParameter{carVar}(25 / 180 * pi, " max steering angle", "rad")
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
        lessContraint(steeringAngle.value / maxAngle.value, 1.0, model)
        greaterContraint(steeringAngle.value / maxAngle.value, -1.0, model)
    end

    function pivot2CoG(forces::Vector{carVar})
        #returns moments on cog from wheel forces
        torque.value = cross(position.value, forces)
    end

    function setPivotForce(forcesIn::Vector{carVar})
        forces.value = rotZ(steeringAngle.value) * forcesIn
    end

    function setPivotVelocity(angularVelocity::Vector{carVar}, CoGvelocity::Vector{carVar})
        velocityPivot.value = CoGvelocity + cross(angularVelocity, position.value) #CoG2wheelAssembly(velocity.value,angularVelocity)
    end

    function setTireSpeeds()
        velocityTire.value = rotZinv(steeringAngle.value) * velocityPivot.value
    end

    function setVelocity(angularVelocity::Vector{carVar}, velocity::Vector{carVar})
        setPivotVelocity(angularVelocity, velocity)
        setTireSpeeds()
    end

    function getTorque(forcesIn::Vector{carVar})
        setPivotForce(forcesIn)
        pivot2CoG(forces.value)
    end

    function setupObservables(ax::Axis, tire::Tire)
        return tire.setupObservables(ax)
    end

    function updateObservables(obs::NamedTuple, tire::Tire, x::Float64, y::Float64, ψ::Float64, Fz_static::Float64)
        R = _rotmat2d(ψ)
        global_pos = R * [position.value[1], position.value[2]] + [x, y]
        total_angle = ψ + steeringAngle.value
        tire.updateObservables(obs, global_pos[1], global_pos[2], total_angle, Fz_static)
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
        getTorque,
        setupObservables,
        updateObservables,
    )
    return testWheelAssembly
end

function createBusWheelAssembly(position::Vector{carVar})
    steeringAngle = carParameter{carVar}(0.0, "steering angle", "rad")
    maxAngle = carParameter{carVar}(45 / 180 * pi, " max steering angle", "rad")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Pivot forces", "N N N")
    torque = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Generated torque on CoG", "N N N")
    velocityPivot = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Velocity at pivot", "m/s")
    velocityTire = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Velocity in tire frame", "m/s")
    position = carParameter{Vector{carVar}}(position, "Position from CoG", "m m m")
    function rotZ(steering::carVar)
        out = [
            cos(steering) -sin(steering) 0;
            sin(steering) cos(steering) 0;
            0 0 1
        ]
        return out
    end

    function rotZinv(steering::carVar)
        out = [
            cos(steering) sin(steering) 0;
            -sin(steering) cos(steering) 0;
            0 0 1
        ]
        return out
    end
    function constraints(model=nothing)
        lessContraint(steeringAngle.value / maxAngle.value, 1.0, model)
        greaterContraint(steeringAngle.value / maxAngle.value, -1.0, model)
    end

    function pivot2CoG(forces::Vector{carVar})
        torque.value = cross(position.value, forces)
    end

    function setPivotForce(forcesIn::Vector{carVar})
        forces.value = rotZ(steeringAngle.value) * forcesIn
    end

    function setPivotVelocity(angularVelocity::Vector{carVar}, CoGvelocity::Vector{carVar})
        velocityPivot.value = CoGvelocity + cross(angularVelocity, position.value)
    end

    function setTireSpeeds()
        velocityTire.value = rotZinv(steeringAngle.value) * velocityPivot.value
    end

    function setVelocity(angularVelocity::Vector{carVar}, velocity::Vector{carVar})
        setPivotVelocity(angularVelocity, velocity)
        setTireSpeeds()
    end

    function getTorque(forcesIn::Vector{carVar})
        setPivotForce(forcesIn)
        pivot2CoG(forces.value)
    end

    function setupObservables(ax::Axis, tire::Tire)
        return tire.setupObservables(ax)
    end

    function updateObservables(obs::NamedTuple, tire::Tire, x::Float64, y::Float64, ψ::Float64, Fz_static::Float64)
        R = _rotmat2d(ψ)
        global_pos = R * [position.value[1], position.value[2]] + [x, y]
        total_angle = ψ + steeringAngle.value
        tire.updateObservables(obs, global_pos[1], global_pos[2], total_angle, Fz_static)
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
        getTorque,
        setupObservables,
        updateObservables,
    )
    return testWheelAssembly
end

"""
    draw!(ax, wa::WheelAssembly, tire::Tire, x, y, ψ)

Draw the wheel assembly by drawing its tire at the correct global position.
Uses the assembly's position offset and steering angle, plus the car's CoG (x,y) and heading ψ.
"""
function draw!(ax, wa::WheelAssembly, tire::Tire, x, y, ψ)
    R = _rotmat2d(ψ)
    global_pos = R * [wa.position.value[1], wa.position.value[2]] + [x, y]
    total_angle = ψ + wa.steeringAngle.value
    draw!(ax, tire, global_pos[1], global_pos[2], total_angle)
end

