using SLapSim
function createBus(track::Union{Track,Nothing}=nothing)

    if isnothing(track)
        widthR = 3.0
        widthL = 3.0
    else
        widthL = maximum(track.widthL)
        widthR = maximum(track.widthR)
    end

    velocity        = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s", :static, [2.0, 30.0])
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s", :static, [-5.0, 5.0])

    mass            = carParameter{carVar}(12000.0, "Mass", "kg")
    motorForce      = carParameter{carVar}(3000.0, "motorForce", "N")
    lateralForce    = carParameter{carVar}(0.0, "lateral Force", "N")
    lateralTransfer = carParameter{carVar}(0.0, "lateral load transfer", "N")
    brakeBias       = carParameter{carVar}(0.6, "brake bias front", "-", :sensitivity)
    CL              = carParameter{carVar}(0.0, "Lift Coefficient", "-")
    CD              = carParameter{carVar}(6.0, "Drag Coefficient", "-")
    powerLimit      = carParameter{carVar}(300000.0, "PowerLimit", "W")
    psi             = carParameter{carVar}(0.0, "heading", "rad", :static, [-4pi, 4pi])
    n               = carParameter{carVar}(0.0, "Distance from centerline", "m", :static, [-widthR+1.5, widthL-1.5])
    brakeCommand    = carParameter{carVar}(0.0, "brake command", "N", :control, [-500000.0, 0.0])
    nControls       = carParameter{carVar}(3.0, "number of controlled parameters", "-")
    inertia         = carParameter{carVar}(30000.0, "Inertia", "kg*m^2")
    nStates         = carParameter{carVar}(6.0, "number of car states", "-")
    s               = carParameter{carVar}(1.0, "longitudinal position on track", "-", :static, [0.0, 200.0])

    # rear axle spacing: 1.3m between the two rear axles
    rearAxleSpacing = 1.3

    # 6 gearboxes from formula, override ratio
    gearboxes = [createCTU25gearbox() for _ in 1:6]
    for gb in gearboxes
        gb.ratio.value = 6.0
    end

    # 6 motors from formula with bus torque limit
    motors = [createFischerMotor(1000.0) for _ in 1:6]
    for m in motors
        m.mass.value = 120.0
    end

    # 6 tires from formula, override dimensions
    tires = [createR20_pacejka(motors[i], gearboxes[i]) for i in 1:6]
    for t in tires
        t.radius.value = 0.5
        t.width.value = 0.315
        t.inertia.value = 5.0
        t.mass.value = 50.0
        t.maxSLipAngle.value = 8/180*pi
        t.scalingForce.value = motors[1].torqueSpeedFunction(0.0) * gearboxes[1].ratio.value / 0.5
    end

    # accumulator from formula, override capacity
    accu = createAccumulator()
    accu.capacity.value = 200.0
    accu.maxPower.value = 300.0
    accu.minPower.value = 300.0
    accu.voltage.value = 650.0
    accu.mass.value = 1500.0

    drivetrain = Drivetrain(motors, gearboxes, tires, accu)

    aero = createBasicAero()
    aero.CL.value = -0.2
    aero.CD.value = 6.0

    suspension = createBusSuspension()

    chassis = createCTU25chassis()
    chassis.mass.value = 12000.0
    chassis.wheelbase.value = 6.0
    chassis.track.value = 2.1
    chassis.CoG_X_pos.value = 0.45

    suspension.setInput(chassis)

    # wheel positions relative to CoG
    cogOffsetX = (chassis.CoG_X_pos.value - 0.5) * chassis.wheelbase.value
    cogOffsetY = (chassis.CoG_Y_pos.value - 0.5) * chassis.track.value
    rearCenter = -chassis.wheelbase.value / 2 - cogOffsetX
    wheelAssemblyFL  = createBusWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX,  chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyFR  = createBusWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRL1 = createBasicWheelAssembly(Vector{carVar}([rearCenter + rearAxleSpacing / 2,  chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRR1 = createBasicWheelAssembly(Vector{carVar}([rearCenter + rearAxleSpacing / 2, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRL2 = createBasicWheelAssembly(Vector{carVar}([rearCenter - rearAxleSpacing / 2,  chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRR2 = createBasicWheelAssembly(Vector{carVar}([rearCenter - rearAxleSpacing / 2, -chassis.track.value / 2 - cogOffsetY, 0]))

    wheelAssemblies = [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL1, wheelAssemblyRR1, wheelAssemblyRL2, wheelAssemblyRR2]

    state_descriptor = VarEntry[
        VarEntry("vx",    [velocity => 1]),
        VarEntry("vy",    [velocity => 2],          -5.0,  5.0),
        VarEntry("psi",   [psi => 0]),
        VarEntry("omega", [angularVelocity => 3]),
        VarEntry("n",     [n => 0]),
        VarEntry("t",     [s => 0]),
    ]

    control_descriptor = VarEntry[
        VarEntry("torque", [
            #drivetrain.motors[1].torque => 0,
            #drivetrain.motors[2].torque => 0,
            drivetrain.motors[3].torque => 0,
            drivetrain.motors[4].torque => 0,
            #drivetrain.motors[5].torque => 0,
            #drivetrain.motors[6].torque => 0,
        ]),
        VarEntry("steering_front", [wheelAssemblies[1].steeringAngle => 0, wheelAssemblies[2].steeringAngle => 0]),
        VarEntry("steering_rear",  [wheelAssemblies[5].steeringAngle => 0, wheelAssemblies[6].steeringAngle => 0]),
        VarEntry("brake", [brakeCommand => 0]),
    ]

    nControls.value = Float64(length(control_descriptor))
    nStates.value   = Float64(length(state_descriptor))

    function controlMapping(controls::AbstractVector)
        apply_mapping!(control_descriptor, controls)
        frontRatio = chassis.CoG_X_pos.value
        rearRatio  = 1 - chassis.CoG_X_pos.value
        leftRatio  = 1 - chassis.CoG_Y_pos.value
        rightRatio = chassis.CoG_Y_pos.value

        drivetrain.tires[1].brakingForce.value = brakeCommand.value * frontRatio * leftRatio           # FL
        drivetrain.tires[2].brakingForce.value = brakeCommand.value * frontRatio * rightRatio          # FR
        drivetrain.tires[3].brakingForce.value = brakeCommand.value * rearRatio  * leftRatio  / 2      # RL1
        drivetrain.tires[4].brakingForce.value = brakeCommand.value * rearRatio  * rightRatio / 2      # RR1
        drivetrain.tires[5].brakingForce.value = brakeCommand.value * rearRatio  * leftRatio  / 2      # RL2
        drivetrain.tires[6].brakingForce.value = brakeCommand.value * rearRatio  * rightRatio / 2      # RR2
    end

    function stateMapping(states::AbstractVector)
        apply_mapping!(state_descriptor, states)
    end

    function anyTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)

        for wa in wheelAssemblies
            wa.setVelocity(angularVelocity.value, velocity.value)
        end

        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].setVelocity(wheelAssemblies[i].velocityTire.value)
        end

        for i in eachindex(drivetrain.gearboxes)
            drivetrain.gearboxes[i].setTorque(drivetrain.motors[i].torque.value)
            drivetrain.gearboxes[i].compute()
        end

        aeroForces = aero.compute(velocity.value[1], isnothing(track) ? RHO_SEA_LEVEL : track.rho[1])

        forces = suspension.calculate(aeroForces.downforce, aero.CoP.value)
        for i in eachindex(forces)
            drivetrain.tires[i].forces.value[3] = forces[i]
        end

        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].compute(drivetrain.gearboxes[i].torqueOut.value, optiModel)
        end

        ######################################################CONSTRAINTS###########################
        for i in eachindex(drivetrain.motors)
            #drivetrain.motors[i].constraints(drivetrain.motors[i].torque.value, optiModel)
        end
        #wheelAssemblies[1].constraints(optiModel)
        chassis.hitbox(n.value, track, optiModel)
        for tire in drivetrain.tires
            tire.tireConstraints(optiModel)
        end

        for i in eachindex(wheelAssemblies)
            wheelAssemblies[i].getTorque(drivetrain.tires[i].forces.value)
        end

        cogForce  = zero(wheelAssemblies[1].forces.value)
        cogMoment = zero(wheelAssemblies[1].forces.value)
        for i in eachindex(wheelAssemblies)
            cogForce  = cogForce  .+ wheelAssemblies[i].forces.value
            cogMoment = cogMoment  + wheelAssemblies[i].torque.value
        end

        cogForce = cogForce .+ [aeroForces.drag, 0.0, 0.0]
        dv = cogForce ./ mass.value - angularVelocity.value × velocity.value
        dangularVelocity = cogMoment ./ inertia.value
        dx = [dv[1], dv[2], angularVelocity.value[3], dangularVelocity[3]]
        return dx
    end

    p = CarParameters(
        mass,
        inertia,
        motorForce,
        CL,
        CD,
        velocity,
        angularVelocity,
        psi,
        n,
        powerLimit,
        lateralForce,
        lateralTransfer,
        brakeBias,
        nControls,
        nStates,
        s,
        state_descriptor,
        control_descriptor,
    )

    car = Car(
        anyTrack,
        p,
        controlMapping,
        stateMapping,
        f -> (0.0),
        drivetrain,
        aero,
        suspension,
        chassis,
        wheelAssemblies,
        state_descriptor,
        control_descriptor,
    )
    return car
end
