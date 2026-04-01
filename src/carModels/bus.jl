using Revise
using SLapSim
function createBus()
    velocity = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s")
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s")

    mass = carParameter{carVar}(12000.0, "Mass", "kg")
    motorForce = carParameter{carVar}(3000.0, "motorForce", "N")
    lateralForce = carParameter{carVar}(0.0, "lateral Force", "N")
    CL = carParameter{carVar}(0.0, "Lift Coefficient", "-")
    CD = carParameter{carVar}(6.0, "Drag Coefficient", "-")
    powerLimit = carParameter{carVar}(300000.0, "PowerLimit", "W")
    psi = carParameter{carVar}(0.0, "heading", "rad")
    n = carParameter{carVar}(0.0, "Distance from centerline", "m")
    nControls = carParameter{carVar}(2.0, "number of controlled parameters", "-")
    inertia = carParameter{carVar}(30000.0, "Inertia", "kg*m^2")
    nStates = carParameter{carVar}(6.0, "number of car states", "-")
    s = carParameter{carVar}(1.0, "longitudinal position on track", "-")

    # rear axle spacing: 1.3m between the two rear axles
    rearAxleSpacing = 1.3

    # 6 gearboxes: FL, FR, RL1, RR1, RL2, RR2
    gearboxFL  = createBusGearbox()
    gearboxFR  = createBusGearbox()
    gearboxRL1 = createBusGearbox()
    gearboxRR1 = createBusGearbox()
    gearboxRL2 = createBusGearbox()
    gearboxRR2 = createBusGearbox()

    # 6 motors
    motorFL  = createBusMotor()
    motorFR  = createBusMotor()
    motorRL1 = createBusMotor()
    motorRR1 = createBusMotor()
    motorRL2 = createBusMotor()
    motorRR2 = createBusMotor()

    # 6 tires
    tireFL  = createBusTire(motorFL, gearboxFL)
    tireFR  = createBusTire(motorFR, gearboxFR)
    tireRL1 = createBusTire(motorRL1, gearboxRL1)
    tireRR1 = createBusTire(motorRR1, gearboxRR1)
    tireRL2 = createBusTire(motorRL2, gearboxRL2)
    tireRR2 = createBusTire(motorRR2, gearboxRR2)

    drivetrain = Drivetrain(
        [motorFL, motorFR, motorRL1, motorRR1, motorRL2, motorRR2],
        [gearboxFL, gearboxFR, gearboxRL1, gearboxRR1, gearboxRL2, gearboxRR2],
        [tireFL, tireFR, tireRL1, tireRR1, tireRL2, tireRR2],
        createBusAccumulator())

    aero = createBusAero()
    suspension = createBusSuspension()

    chassis = createBusChassis()
    suspension.setInput(chassis)

    # wheel positions: front axle, rear axle 1, rear axle 2
    rearCenter = -chassis.wheelbase.value / 2
    wheelAssemblyFL  = createBusWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2, chassis.track.value / 2, 0]))
    wheelAssemblyFR  = createBusWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2, -chassis.track.value / 2, 0]))
    wheelAssemblyRL1 = createBasicWheelAssembly(Vector{carVar}([rearCenter + rearAxleSpacing / 2, chassis.track.value / 2, 0]))
    wheelAssemblyRR1 = createBasicWheelAssembly(Vector{carVar}([rearCenter + rearAxleSpacing / 2, -chassis.track.value / 2, 0]))
    wheelAssemblyRL2 = createBasicWheelAssembly(Vector{carVar}([rearCenter - rearAxleSpacing / 2, chassis.track.value / 2, 0]))
    wheelAssemblyRR2 = createBasicWheelAssembly(Vector{carVar}([rearCenter - rearAxleSpacing / 2, -chassis.track.value / 2, 0]))

    wheelAssemblies = [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL1, wheelAssemblyRR1, wheelAssemblyRL2, wheelAssemblyRR2]

    function controlMapping(controls::AbstractVector)
        # front motors off
        drivetrain.motors[1].torque.value = 0.0
        drivetrain.motors[2].torque.value = 0.0
        # rear motors - all 4 rear wheels get same torque (controls[1])
        drivetrain.motors[3].torque.value = controls[1]
        drivetrain.motors[4].torque.value = controls[1]
        #drivetrain.motors[5].torque.value = controls[1]
        #drivetrain.motors[6].torque.value = controls[1]
        # front steering only
        wheelAssemblies[1].steeringAngle.value = controls[2]
        wheelAssemblies[2].steeringAngle.value = controls[2]
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

        forces = suspension.calculate()
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
