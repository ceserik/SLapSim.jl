using SLapSim
function createBus()
    velocity = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s")
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s")

    mass = carParameter{carVar}(12000.0, "Mass", "kg")
    motorForce = carParameter{carVar}(3000.0, "motorForce", "N")
    lateralForce = carParameter{carVar}(0.0, "lateral Force", "N")
    lateralTransfer = carParameter{carVar}(0.0, "lateral load transfer", "N")
    brakeBias = carParameter{carVar}(0.6, "brake bias front", "-", :tunable)
    CL = carParameter{carVar}(0.0, "Lift Coefficient", "-")
    CD = carParameter{carVar}(6.0, "Drag Coefficient", "-")
    powerLimit = carParameter{carVar}(300000.0, "PowerLimit", "W")
    psi = carParameter{carVar}(0.0, "heading", "rad")
    n = carParameter{carVar}(0.0, "Distance from centerline", "m")
    nControls = carParameter{carVar}(3.0, "number of controlled parameters", "-")
    inertia = carParameter{carVar}(30000.0, "Inertia", "kg*m^2")
    nStates = carParameter{carVar}(6.0, "number of car states", "-")
    s = carParameter{carVar}(1.0, "longitudinal position on track", "-")

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
    tires = [createR20lin(motors[i], gearboxes[i]) for i in 1:6]
    for t in tires
        t.radius.value = 0.5
        t.width.value = 0.315
        t.inertia.value = 5.0
        t.mass.value = 50.0
        t.maxSLipAngle.value = 8/180*pi
        t.scalingForce.value = motors[1].torqueSpeedFunction(0.0) * gearboxes[1].ratio.value / 0.5
    end

    # accumulator from formula, override capacity
    accu = createPepikCTU25()
    accu.capacity.value = 200.0
    accu.maxPower.value = 300.0
    accu.minPower.value = 300.0
    accu.voltage.value = 650.0
    accu.mass.value = 1500.0

    drivetrain = Drivetrain(motors, gearboxes, tires, accu)

    # aero from formula, override for bus (no downforce, high drag)
    aero = createBasicAero()
    aero.CL.value = 0.0
    aero.CD.value = 6.0

    # suspension - use simple suspension but returns 6 forces for double rear axle
    suspension = createBusSuspension()

    # chassis from formula, override dimensions
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
    wheelAssemblyFL  = createBusWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyFR  = createBusWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRL1 = createBasicWheelAssembly(Vector{carVar}([rearCenter + rearAxleSpacing / 2, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRR1 = createBasicWheelAssembly(Vector{carVar}([rearCenter + rearAxleSpacing / 2, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRL2 = createBasicWheelAssembly(Vector{carVar}([rearCenter - rearAxleSpacing / 2, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRR2 = createBasicWheelAssembly(Vector{carVar}([rearCenter - rearAxleSpacing / 2, -chassis.track.value / 2 - cogOffsetY, 0]))

    wheelAssemblies = [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL1, wheelAssemblyRR1, wheelAssemblyRL2, wheelAssemblyRR2]

    function controlMapping(controls::AbstractVector)
        # front motors off
        drivetrain.motors[1].torque.value = controls[1]
        drivetrain.motors[2].torque.value = controls[1]
        # rear motors - all 4 rear wheels get same torque (controls[1])
        drivetrain.motors[3].torque.value = controls[1]
        drivetrain.motors[4].torque.value = controls[1]
        drivetrain.motors[5].torque.value = controls[1]
        drivetrain.motors[6].torque.value = controls[1]
        # front steering only
        wheelAssemblies[1].steeringAngle.value = controls[2]
        wheelAssemblies[2].steeringAngle.value = controls[2]

        wheelAssemblies[5].steeringAngle.value = controls[3]
        wheelAssemblies[6].steeringAngle.value = controls[3]
    end

    function stateMapping(states::AbstractVector)
        velocity.value .= [states[1], states[2], 0.0]
        psi.value = states[3]
        angularVelocity.value .= [0.0, 0.0, states[4]]
        n.value = states[5]
        s.value = states[6]
        return car
    end

    function anyTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)

        for wa in car.wheelAssemblies
            wa.setVelocity(angularVelocity.value, velocity.value)
        end

        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].setVelocity(wheelAssemblies[i].velocityTire.value)
        end

        for i in eachindex(car.drivetrain.gearboxes)
            car.drivetrain.gearboxes[i].setTorque(car.drivetrain.motors[i].torque.value)
            car.drivetrain.gearboxes[i].compute()
        end

        aeroForces = aero.compute(velocity.value[1], isnothing(track) ? RHO_SEA_LEVEL : track.rho[1])

        forces = suspension.calculate(aeroForces.downforce, aero.CoP.value)
        for i in eachindex(forces)
            drivetrain.tires[i].forces.value[3] = forces[i]
        end

        for i in eachindex(car.drivetrain.tires)
            car.drivetrain.tires[i].compute(car.drivetrain.gearboxes[i].torqueOut.value, optiModel)
        end

######################################################CONSTRAINTS###########################
        # constrain all motors
        for i in eachindex(car.drivetrain.motors)
            car.drivetrain.motors[i].constraints(drivetrain.motors[i].torque.value, optiModel)
        end
        # front steering constraint
        car.wheelAssemblies[1].constraints(optiModel)
        # hitbox
        car.chassis.hitbox(car.carParameters.n.value, track, optiModel)
        for tire in car.drivetrain.tires
            tire.tireConstraints(optiModel)
        end

        for i in eachindex(car.wheelAssemblies)
            car.wheelAssemblies[i].getTorque(car.drivetrain.tires[i].forces.value)
        end

        cogForce = zero(car.wheelAssemblies[1].forces.value)
        cogMoment = zero(car.wheelAssemblies[1].forces.value)
        for i in eachindex(car.wheelAssemblies)
            cogForce = cogForce .+ car.wheelAssemblies[i].forces.value
            cogMoment = cogMoment + car.wheelAssemblies[i].torque.value
        end

        cogForce = cogForce .+ [aeroForces.drag, 0.0, 0.0]
        dv = cogForce / car.carParameters.mass.value - angularVelocity.value × velocity.value
        dangularVelocity = cogMoment / car.carParameters.inertia.value
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
        s
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
        [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL1, wheelAssemblyRR1, wheelAssemblyRL2, wheelAssemblyRR2],
    )
    return car
end
