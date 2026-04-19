function createTwintrack(pacejka::Bool=true,track::Union{Track,Nothing} = nothing)

    velocity         = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s", :static, [3.0, 60.0])
    angularVelocity  = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s", :static, [-10.0, 10.0])

    mass             = carParameter{carVar}(280.0, "Mass", "kg", :tunable, [200.0, 320.0])
    motorForce       = carParameter{carVar}(1000.0, "motorForce", "N")
    lateralForce     = carParameter{carVar}(0.0, "lateral Force", "N")
    lateralTransfer  = carParameter{carVar}(0.0, "lateral load transfer", "N")
    brakeBias        = carParameter{carVar}(0.6, "brake bias front", "-", :tunable)
    CL               = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD               = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    powerLimit       = carParameter{carVar}(80000.0, "PowerLimit", "W")
    psi              = carParameter{carVar}(0.0, "heading", "rad", :static, [-4pi, 4pi])
    n                = carParameter{carVar}(0.0, "Distance from centerline", "m", :static)
    brakeCommand     = carParameter{carVar}(0.0, "brake command", "N", :control, [-20000.0, 0.0])
    nControls        = carParameter{carVar}(0.0, "number of controlled parameters", "-")
    inertia          = carParameter{carVar}(100.0, "Inertia", "kg*m^2", :tunable)
    nStates          = carParameter{carVar}(0.0, "number of car states", "-")
    s                = carParameter{carVar}(1.0, "longitudinal position on track", "-", :static, [0.0, 200.0])

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
    accumulator = createAccumulator()
    drivetrain = Drivetrain(
        [motorFL, motorFR, motorRL, motorRR],
        [gearboxFL, gearboxFR, gearboxRL, gearboxRR],
        [tireFL, tireFR, tireRL, tireRR],
        accumulator)

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

    state_descriptor = VarEntry[
        VarEntry("vx",    [velocity => 1],          :state),
        VarEntry("vy",    [velocity => 2], -1.5, 1.5, :state),
        VarEntry("psi",   [psi => 0],               :state),
        VarEntry("omega", [angularVelocity => 3],   :state),
        VarEntry("n",     [n => 0],                 :state),
        VarEntry("t",     [s => 0],                 :state),
    ]

    control_descriptor = VarEntry[
        VarEntry("torque_rear",  [drivetrain.motors[3].torque => 0, drivetrain.motors[4].torque => 0], :control),
        VarEntry("steering",    [wheelAssemblies[1].steeringAngle => 0, wheelAssemblies[2].steeringAngle => 0], :control),
        VarEntry("torque_front", [drivetrain.motors[1].torque => 0, drivetrain.motors[2].torque => 0], :control),
        VarEntry("brake", [brakeCommand => 0], :control),
    ]

    nControls.value = Float64(length(control_descriptor))
    nStates.value   = Float64(length(state_descriptor))

    function controlMapping(controls::AbstractVector)
        apply_mapping!(control_descriptor, controls)
        apply_bounds!(control_descriptor)
        frontBrake = brakeCommand.value * brakeBias.value
        rearBrake = brakeCommand.value * (1 - brakeBias.value)
        drivetrain.tires[1].brakingForce.value = frontBrake
        drivetrain.tires[2].brakingForce.value = frontBrake
        drivetrain.tires[3].brakingForce.value = rearBrake
        drivetrain.tires[4].brakingForce.value = rearBrake
    end

    function stateMapping(states::AbstractVector)
        apply_mapping!(state_descriptor, states)
    end

    function anyTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)
        # Transformation of velocities from cog to wheels and steering
        if !isnothing(optiModel)
            # Prevent simultaneous acceleration and braking
            T_max = drivetrain.motors[3].torque.limits[2]  # 29 Nm
            F_max = -brakeCommand.limits[1]                # 20000 N
            # complementarity brake and motors cannot go againts each other
            @constraint(optiModel, (drivetrain.motors[3].torque.value / T_max) * (-brakeCommand.value / F_max) <= 0.001)
            @constraint(optiModel, (drivetrain.motors[1].torque.value / T_max) * (-brakeCommand.value / F_max) <= 0.001)
        end
        for wa in wheelAssemblies
            wa.setVelocity(angularVelocity.value, velocity.value)
        end

        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].setVelocity(wheelAssemblies[i].velocityTire.value)
        end

        #gearing of forces from motor to tire
        for i in eachindex(drivetrain.gearboxes)
            drivetrain.gearboxes[i].setTorque(drivetrain.motors[i].torque.value)
            drivetrain.gearboxes[i].compute()
        end

        # gearing of tire velocity to motor velocity
        for i in eachindex(drivetrain.motors)
            drivetrain.motors[i].setVelocity(drivetrain.tires[i].angularFrequency.value * drivetrain.gearboxes[i].ratio.value)
        end

        drivetrain.accumulators.compute(drivetrain.motors)
        #compute aero
        aeroForces = aero.compute(velocity.value[1], isnothing(track) ? RHO_SEA_LEVEL : track.rho[1])

        #compute suspension forces
        forces = suspension.calculate(aeroForces.downforce, aero.CoP.value)
        for i in eachindex(forces)
            drivetrain.tires[i].forces.value[3] = forces[i]
        end

        #tire forces
        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].compute(drivetrain.gearboxes[i].torqueOut.value, optiModel)
        end

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
        s,
        state_descriptor,
        control_descriptor
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
        state_descriptor,
        control_descriptor,
    )
    return afto
end


function formulaE2026(track::Union{Track,Nothing}=nothing)


    velocity         = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s", :static, [2.0, 60.0])
    angularVelocity  = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s", :static, [-10.0, 10.0])
    mass             = carParameter{carVar}(1200.0, "Mass", "kg", :tunable)
    motorForce       = carParameter{carVar}(7100.0, "motorForce", "N")
    lateralForce     = carParameter{carVar}(0.0, "lateral Force", "N")
    brakeBias        = carParameter{carVar}(0.7, "brake bias front", "-", :tunable)
    CL               = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD               = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    powerLimit       = carParameter{carVar}(270000.0, "PowerLimit", "W")
    psi              = carParameter{carVar}(0.0, "heading", "rad", :static, [-4pi, 4pi])
    n                = carParameter{carVar}(0.0, "Distance from centerline", "m", :static)
    inertia          = carParameter{carVar}(1260.0, "Inertia", "kg*m^2", :tunable)
    s                = carParameter{carVar}(1.0, "longitudinal position on track", "-", :static, [0.0, 200.0])
    brakeCommand     = carParameter{carVar}(0.0, "brake command", "N", :control, [-20000.0, 0])
    nControls        = carParameter{carVar}(0.0, "number of controlled parameters", "-")
    nStates          = carParameter{carVar}(0.0, "number of car states", "-")

    gearbox = createCTU25gearbox(1.0)
    motor   = createFischerMotor(7100*0.33,0.0)#tire radius to get moment

    tireFL = create_FormulaE_pacejka(motor, gearbox)
    tireFR = create_FormulaE_pacejka(motor, gearbox)
    tireRL = create_FormulaE_pacejka(motor, gearbox)
    tireRR = create_FormulaE_pacejka(motor, gearbox)

    

    drivetrain = Drivetrain(
        [motor],
        [gearbox],
        [tireFL, tireFR, tireRL, tireRR],
        createAccumulator())

    aero = createBasicAero(
        CL_a  = -2.7,
        CD_a  = 1.4,
        CoP_a = 0.44
        )

    chassis = createCTU25chassis(
        mass_p      = 1200.0,
        wheelbase_p = 2.9,
        track_p     = 1.5, CoGx=0.482,
        CogY        = 0.5,
        width_p     = 2.0
        )

    suspension = createQuasi_steady_Suspension()

    cogOffsetX = (chassis.CoG_X_pos.value - 0.5) * chassis.wheelbase.value
    cogOffsetY = (chassis.CoG_Y_pos.value - 0.5) * chassis.track.value
    wheelAssemblyFL = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyFR = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRL = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRR = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblies = [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL, wheelAssemblyRR]


    state_descriptor = VarEntry[
        VarEntry("vx",    [velocity => 1],          :state),
        VarEntry("vy",    [velocity => 2], -5.0, 5.0, :state),
        VarEntry("psi",   [psi => 0],               :state),
        VarEntry("omega", [angularVelocity => 3],   :state),
        VarEntry("n",     [n => 0],                 :state),
        VarEntry("t",     [s => 0],                 :state),
    ]

    control_descriptor = VarEntry[
        VarEntry("torque",           [drivetrain.motors[1].torque => 0],                                              :control),
        VarEntry("steering",         [wheelAssemblies[1].steeringAngle => 0, wheelAssemblies[2].steeringAngle => 0], :control),
        VarEntry("lateral transfer", [suspension.lateralTransfer => 0],                                              :control),
        VarEntry("brake",            [brakeCommand => 0],                                                            :control),
    ]

    nControls.value = Float64(length(control_descriptor))
    nStates.value   = Float64(length(state_descriptor))

    function controlMapping(controls::AbstractVector)
        apply_mapping!(control_descriptor, controls)

        frontBrake = brakeCommand.value * brakeBias.value
        rearBrake = brakeCommand.value * (1 - brakeBias.value)
        drivetrain.tires[1].brakingForce.value = frontBrake
        drivetrain.tires[2].brakingForce.value = frontBrake
        drivetrain.tires[3].brakingForce.value = rearBrake
        drivetrain.tires[4].brakingForce.value = rearBrake
    end

    function stateMapping(states::AbstractVector)
        apply_mapping!(state_descriptor, states)
    end

    function anyTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)
        if !isnothing(optiModel)
            T_max = drivetrain.motors[1].torque.limits[2]
            F_max = -brakeCommand.limits[1]
            scale = T_max * F_max
            @constraint(optiModel, drivetrain.motors[1].torque.value * brakeCommand.value / scale >= -0.0001)
        end
        # Transformation of velocities from cog to wheels and steering
        for wa in wheelAssemblies
            wa.setVelocity(angularVelocity.value, velocity.value)
        end
        # set Tire velocities
        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].setVelocity(wheelAssemblies[i].velocityTire.value)
        end

        # gearing of forces from motor to tire
        drivetrain.gearboxes[1].setTorque(drivetrain.motors[1].torque.value)
        drivetrain.gearboxes[1].compute()


        # Acceleration approximation
        aeroForces = aero.compute(velocity.value[1], isnothing(track) ? RHO_SEA_LEVEL : track.rho[1])
        ax = (drivetrain.gearboxes[1].torqueOut.value / drivetrain.tires[1].radius.value +
        sum(t.brakingForce.value for t in drivetrain.tires) +
        aeroForces[2] + 
        drivetrain.tires[1].rollingResistance.value * chassis.mass.value * 9.81)/chassis.mass.value


        ## suspension calc
        suspension.setInput(chassis,wheelAssemblies)
        forces = suspension.calculate(ax,aeroForces[1],aero.CoP.value,optiModel)
        for i in eachindex(forces)
            drivetrain.tires[i].forces.value[3] = forces[i]
        end

        #tire forces
        torqueOut = drivetrain.gearboxes[1].torqueOut.value
        #for i in eachindex(drivetrain.tires)


        drivetrain.tires[1].compute(0.0, optiModel)
        drivetrain.tires[2].compute(0.0, optiModel)
        drivetrain.tires[3].compute(torqueOut/2, optiModel)
        drivetrain.tires[4].compute(torqueOut/2, optiModel)
        #end

        for tire in drivetrain.tires
            tire.tireConstraints(optiModel)
        end

        for i in eachindex(wheelAssemblies)
            wheelAssemblies[i].getTorque(drivetrain.tires[i].forces.value)
        end


        if !isnothing(optiModel)
            @constraint(optiModel, suspension.lateralTransfer.value == (chassis.CoG_Z_pos.value / chassis.track.value) * sum(wa.forces.value[2] for wa in wheelAssemblies))
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
        suspension.lateralTransfer,
        brakeBias,
        nControls,
        nStates,
        s,
        state_descriptor,
        control_descriptor
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
        nothing,
        nothing,
    )
    return afto

end