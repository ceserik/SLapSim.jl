function createTwintrack(pacejka::Bool=true)
    velocity = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s")
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s")

    mass = carParameter{carVar}(280.0, "Mass", "kg", :tunable)
    motorForce = carParameter{carVar}(1000.0, "motorForce", "N")
    lateralForce = carParameter{carVar}(0.0, "lateral Force", "N")
    lateralTransfer = carParameter{carVar}(0.0, "lateral load transfer", "N")
    brakeBias = carParameter{carVar}(0.6, "brake bias front", "-", :tunable)
    CL = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    powerLimit = carParameter{carVar}(80000.0, "PowerLimit", "W")
    psi = carParameter{carVar}(0.0, "heading", "rad")
    n = carParameter{carVar}(0.0, "Distance from centerline", "m")
    nControls = carParameter{carVar}(3.0, "number of controlled parameters", "-")
    inertia = carParameter{carVar}(100.0, "Inertia", "kg*m^2", :tunable)
    nStates = carParameter{carVar}(6.0, "number of car states", "-")
    s = carParameter{carVar}(1.0, "longitudinal position on track", "-")

    gearboxFL = createCTU25gearbox()
    gearboxFR = createCTU25gearbox()
    gearboxRL = createCTU25gearbox()
    gearboxRR = createCTU25gearbox()

    motorFL = createFischerMotor()
    motorFR = createFischerMotor()
    motorRL = createFischerMotor()
    motorRR = createFischerMotor()

    if pacejka
        tireFL = createR20_pacejka(motorFL, gearboxFL)
        tireFR = createR20_pacejka(motorFR, gearboxFR)
        tireRL = createR20_pacejka(motorRL, gearboxRL)
        tireRR = createR20_pacejka(motorRR, gearboxRR)
    else
        tireFL = createR20lin(motorFL, gearboxFL)
        tireFR = createR20lin(motorFR, gearboxFR)
        tireRL = createR20lin(motorRL, gearboxRL)
        tireRR = createR20lin(motorRR, gearboxRR)

    end
    drivetrain = Drivetrain(
        [motorFL, motorFR, motorRL, motorRR],
        [gearboxFL, gearboxFR, gearboxRL, gearboxRR],
        [tireFL, tireFR, tireRL, tireRR],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createSimpleSuspension()

    chassis = createCTU25chassis()
    suspension.setInput(chassis)
    cogOffsetX = (chassis.CoG_X_pos.value - 0.5) * chassis.wheelbase.value
    cogOffsetY = (chassis.CoG_Y_pos.value - 0.5) * chassis.track.value
    wheelAssemblyFL = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyFR = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRL = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRR = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))

    wheelAssemblies = [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL, wheelAssemblyRR]

    function controlMapping(controls::AbstractVector)
        drivetrain.motors[1].torque.value = controls[3]
        drivetrain.motors[2].torque.value = controls[3]
        drivetrain.motors[3].torque.value = controls[1]
        drivetrain.motors[4].torque.value = controls[1]
        wheelAssemblies[1].steeringAngle.value = controls[2]
        wheelAssemblies[2].steeringAngle.value = controls[2]

    end

    function stateMapping(states::AbstractVector)
        #@infiltrate
        velocity.value .= [states[1], states[2], 0.0]
        psi.value = states[3]
        angularVelocity.value .= [0.0, 0.0, states[4]]
        n.value = states[5]
        s.value = states[6]

    end

    function anyTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)

        #velocity = car.carParameters.velocity.value
        #angularVelocity = car.carParameters.angularVelocity.value
        # Transformation of velocities from cog to wheels and steering
        for wa in wheelAssemblies
            wa.setVelocity(angularVelocity.value, velocity.value)
            wa.constraints(nothing)
        end

        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].setVelocity(wheelAssemblies[i].velocityTire.value)
        end

        #gearing of forces from motor to tire
        for i in eachindex(drivetrain.gearboxes)
            drivetrain.gearboxes[i].setTorque(drivetrain.motors[i].torque.value)
            drivetrain.gearboxes[i].compute()
        end



        aeroForces = aero.compute(velocity.value[1], isnothing(track) ? RHO_SEA_LEVEL : track.rho[1])

        forces = suspension.calculate(aeroForces.downforce, aero.CoP.value)
        for i in eachindex(forces)
            drivetrain.tires[i].forces.value[3] = forces[i]
        end

        #tire forces
        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].compute(drivetrain.gearboxes[i].torqueOut.value, optiModel)
        end

        ######################################################CONSTRAINTS###########################
        # motor torque limit — one per shared control
        drivetrain.motors[1].constraints(drivetrain.motors[1].torque.value, optiModel) # front (controls[3])
        drivetrain.motors[3].constraints(drivetrain.motors[3].torque.value, optiModel) # rear  (controls[1])
        #steering angle
        wheelAssemblies[1].constraints(optiModel)
        #hitbox
        chassis.hitbox(n.value, track, optiModel)
        for tire in drivetrain.tires
            tire.tireConstraints(optiModel)
        end

        for i in eachindex(wheelAssemblies)

            wheelAssemblies[i].getTorque(drivetrain.tires[i].forces.value)
        end


        cogForce = zero(wheelAssemblies[1].forces.value)
        cogMoment = zero(wheelAssemblies[1].forces.value)
        for i in eachindex(wheelAssemblies)
            cogForce = cogForce .+ wheelAssemblies[i].forces.value
            cogMoment = cogMoment + wheelAssemblies[i].torque.value
        end
        cogForce = cogForce .+ [aeroForces.drag, 0.0, 0.0]

        dv = cogForce ./ mass.value - angularVelocity.value × velocity.value #really check what sign should be here !!!! podla mna bednarik skripta fyzika1 Kapitola 8  Neinerciální vztažné soustavy, neboli pro vyjádření časové změny libovolné vektorové veličiny v nečárkované soustavě je možné použít následujícího operátoru:  d·  dt = d′·  dt + ω × · , (8.11)
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
        s
    )

    afto = Car(
        anyTrack,
        p,
        controlMapping,
        stateMapping,
        f -> (0.0),
        drivetrain,
        aero,
        suspension,
        chassis,
        [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL, wheelAssemblyRR],
    )
    return afto
end


function formulaE2026()
    velocity = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s")
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s")
    mass = carParameter{carVar}(1200.0, "Mass", "kg", :tunable)
    motorForce = carParameter{carVar}(1000.0, "motorForce", "N")
    lateralForce = carParameter{carVar}(0.0, "lateral Force", "N")
    lateralTransfer = carParameter{carVar}(0.0, "lateral load transfer", "N")
    brakeBias = carParameter{carVar}(0.6, "brake bias front", "-", :tunable)
    CL = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    powerLimit = carParameter{carVar}(270000.0, "PowerLimit", "W")
    psi = carParameter{carVar}(0.0, "heading", "rad")
    n = carParameter{carVar}(0.0, "Distance from centerline", "m")
    nControls = carParameter{carVar}(4.0, "number of controlled parameters", "-")
    inertia = carParameter{carVar}(1260.0, "Inertia", "kg*m^2", :tunable)
    nStates = carParameter{carVar}(6.0, "number of car states", "-")
    s = carParameter{carVar}(1.0, "longitudinal position on track", "-")



    gearbox = createCTU25gearbox(1.0)
    motor = createFischerMotor()

    tireFL = createR20_pacejka(motor, gearbox)
    tireFR = createR20_pacejka(motor, gearbox)
    tireRL = createR20_pacejka(motor, gearbox)
    tireRR = createR20_pacejka(motor, gearbox)

    drivetrain = Drivetrain(
        [motor],
        [gearbox],
        [tireFL, tireFR, tireRL, tireRR],
        createPepikCTU25())

    aero = createBasicAero(-2.7, 1.4, 0.44)
    chassis = createCTU25chassis(1200.0, 2.9, 1.5, 0.517, 0.5, 2.0) #the model on website has different track for each axle
    suspension = createQuasi_steady_Suspension()

    cogOffsetX = (chassis.CoG_X_pos.value - 0.5) * chassis.wheelbase.value
    cogOffsetY = (chassis.CoG_Y_pos.value - 0.5) * chassis.track.value
    wheelAssemblyFL = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyFR = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRL = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRR = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))

    wheelAssemblies = [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL, wheelAssemblyRR]


    function controlMapping(controls::AbstractVector)
        ## steering
        wheelAssemblies[1].steeringAngle.value = controls[1]
        wheelAssemblies[2].steeringAngle.value = controls[1]
        ## accelerating
        drivetrain.motors[3].torque.value = controls[2]
        drivetrain.motors[4].torque.value = controls[2]
        ## Braking
        frontBrake = controls[3] * brakeBias.value
        rearBrake = controls[3] * (1 - brakeBias.value)
        car.drivetrain.tires[1].brakingForce.value = frontBrake
        car.drivetrain.tires[2].brakingForce.value = frontBrake
        car.drivetrain.tires[3].brakingForce.value = rearBrake
        car.drivetrain.tires[4].brakingForce.value = rearBrake
        # lateral load transfer
        car.suspension.lateralTransfer.value = controls[4]
    end

    function stateMapping(states::AbstractVector)
        velocity.value .= [states[1], states[2], 0.0]
        psi.value = states[3]
        angularVelocity.value .= [0.0, 0.0, states[4]]
        n.value = states[5]
        s.value = states[6]
    end

    function anyTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)

        #velocity = car.carParameters.velocity.value
        #angularVelocity = car.carParameters.angularVelocity.value

        # Transformation of velocities from cog to wheels and steering
        for wa in wheelAssemblies
            wa.setVelocity(angularVelocity.value, velocity.value)
            wa.constraints(optiModel)
        end
        # set Tire velocities
        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].setVelocity(wheelAssemblies[i].velocityTire.value)
        end

        #gearing of forces from motor to tire
        for i in eachindex(drivetrain.gearboxes)
            drivetrain.gearboxes[i].setTorque(drivetrain.motors[i].torque.value)
            drivetrain.gearboxes[i].compute()
        end

        aeroForces = aero.compute(velocity.value[1], isnothing(track) ? RHO_SEA_LEVEL : track.rho[1])
        ax = car.drivetrain.gearboxes[1].torqueOut.value + sum(t.brakingForce.value for t in drivetrain.tires) + aeroForces[2] + drivetrain.tires[1].rollingResistance.value * chassis.mass.value * 9.81

        forces = suspension.calculate(ax)
        for i in eachindex(forces)
            drivetrain.tires[i].forces.value[3] = forces[i]
        end

        #tire forces
        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].compute(drivetrain.gearboxes[i].torqueOut.value, optiModel)
        end

        ######################################################CONSTRAINTS###########################
        # motor torque limit
        drivetrain.motors[1].constraints(drivetrain.motors[1].torque.value, optiModel) # rear (controls[2])
        #steering angle
        wheelAssemblies[1].constraints(optiModel)
        #hitbox
        chassis.hitbox(n.value, track, optiModel)
        for tire in drivetrain.tires
            tire.tireConstraints(optiModel)
        end

        for i in eachindex(wheelAssemblies)

            wheelAssemblies[i].getTorque(drivetrain.tires[i].forces.value)
        end


        cogForce = zero(wheelAssemblies[1].forces.value)
        cogMoment = zero(wheelAssemblies[1].forces.value)
        for i in eachindex(wheelAssemblies)
            cogForce = cogForce .+ wheelAssemblies[i].forces.value
            cogMoment = cogMoment + wheelAssemblies[i].torque.value
        end
        cogForce = cogForce .+ [aeroForces.drag, 0.0, 0.0]

        dv = cogForce ./ mass.value - angularVelocity.value × velocity.value #really check what sign should be here !!!! podla mna bednarik skripta fyzika1 Kapitola 8  Neinerciální vztažné soustavy, neboli pro vyjádření časové změny libovolné vektorové veličiny v nečárkované soustavě je možné použít následujícího operátoru:  d·  dt = d′·  dt + ω × · , (8.11)
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
        s
    )

    afto = Car(
        f -> (0.0),
        p,
        controlMapping,
        stateMapping,
        f -> (0.0),
        drivetrain,
        aero,
        suspension,
        chassis,
        [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL, wheelAssemblyRR],
    )
    return afto

end