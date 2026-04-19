using SLapSim


function createSimplestSingleTrack(track::Union{Track,Nothing}=nothing)

    if isnothing(track)
        widthR = 1.5
        widthL = 1.5
    else
        widthL = maximum(track.widthL)
        widthR = maximum(track.widthR)
    end

    gearboxFront = createCTU25gearbox()
    gearboxRear = createCTU25gearbox()

    motorFront = createFischerMotor()
    motorRear = createFischerMotor()
    
    tireFront = createR20lin(motorFront,gearboxFront)
    tireFront.frictionCoefficient.value = 2.0
    tireRear  = createR20lin(motorRear,gearboxRear)
    tireRear.frictionCoefficient.value = 2.0

    drivetrain = Drivetrain(
        [motorFront, motorRear],
        [gearboxFront, gearboxRear],
        [tireFront, tireRear],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createDummySuspension()
    chassis = createCTU25chassis()
    cogOffsetX = (chassis.CoG_X_pos.value - 0.5) * chassis.wheelbase.value
    wheelAssemblyFront = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, 0, 0]))
    wheelAssemblyRear = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, 0, 0]))

    velocity         = carParameter{Vector{carVar}}([15.0, 0.0, 0.0], "Velocity", "m/s", :static, [3.0, 60.0])
    angularVelocity  = carParameter{Vector{carVar}}([0.0, 0.0, 1.0], "Angular Velocity", "rad/s", :static, [-10.0, 10.0])

    mass = carParameter{carVar}(280.0, "Mass", "kg", :tunable, [200.0, 320.0])
    motorForce = carParameter{carVar}(1000.0, "motorForce", "N")
    lateralForce = carParameter{carVar}(0.0, "lateral Force", "N")
    lateralTransfer = carParameter{carVar}(0.0, "lateral load transfer", "N")
    brakeBias = carParameter{carVar}(0.6, "brake bias front", "-", :tunable)
    CL = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    powerLimit = carParameter{carVar}(80000.0, "PowerLimit", "W")
    psi              = carParameter{carVar}(0.0, "heading", "rad", :static, [-4pi, 4pi])
    n                = carParameter{carVar}(0.0, "Distance from centerline", "m", :static)
    brakeCommand = carParameter{carVar}(0.0, "brake command", "N", :control, [-20000.0, 0.0])
    nControls = carParameter{carVar}(0.0, "number of controlled parameters", "-")
    inertia = carParameter{carVar}(100.0, "Inertia", "kg*m^2", :tunable)
    nStates = carParameter{carVar}(0.0, "number of car states", "-")
    s                = carParameter{carVar}(1.0, "longitudinal position on track", "-", :static, [0.0, 200.0])


    function simplestSingleTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)
        vel = velocity.value
        angVel = angularVelocity.value

        gbFront = drivetrain.gearboxes[1]
        gbRear = drivetrain.gearboxes[2]

        # Transformation of velocities from cog to wheels
        wheelAssemblyFront.setVelocity(angVel,vel)
        wheelAssemblyRear.setVelocity(angVel,vel)
        
        #pass wheel assembly velocities to tires
        drivetrain.tires[1].setVelocity(wheelAssemblyFront.velocityTire.value)
        drivetrain.tires[2].setVelocity(wheelAssemblyRear.velocityTire.value)

        #gearing of forces from motor to tire
        gbFront.setTorque(drivetrain.motors[1].torque.value)
        gbFront.compute()

        gbRear.setTorque(drivetrain.motors[2].torque.value)
        gbRear.compute()

        aeroForces = aero.compute(velocity.value[1], isnothing(track) ? RHO_SEA_LEVEL : track.rho[1])

        # calculate Fz on tires (weight + downforce, split 50/50 for single track)
        drivetrain.tires[1].forces.value[3] = 0.5 * mass.value * 9.81 + aeroForces.downforce * aero.CoP.value
        drivetrain.tires[2].forces.value[3] = 0.5 * mass.value * 9.81 + aeroForces.downforce * (1 - aero.CoP.value)
        
        #tire forces
        drivetrain.tires[1].compute(gbFront.torqueOut.value, optiModel)
        drivetrain.tires[2].compute(gbRear.torqueOut.value, optiModel)

        #if !isnothing(optiModel) ||1==1
        drivetrain.motors[1].constraints(drivetrain.motors[1].torque.value, optiModel)
        drivetrain.motors[2].constraints(drivetrain.motors[2].torque.value, optiModel)
        wheelAssemblyFront.constraints(optiModel)
        wheelAssemblyRear.constraints(optiModel)
        chassis.hitbox(n.value, track, optiModel)
        drivetrain.tires[1].tireConstraints(optiModel)
        drivetrain.tires[2].tireConstraints(optiModel)
        #end

        #propagate forces from tires to wheel assemblies to COG
        wheelAssemblyFront.getTorque(drivetrain.tires[1].forces.value)
        wheelAssemblyRear.getTorque(drivetrain.tires[2].forces.value)

        cogForce = wheelAssemblyFront.forces.value .+ wheelAssemblyRear.forces.value .+ [aeroForces.drag, 0.0, 0.0]
        cogMoment = wheelAssemblyFront.torque.value + wheelAssemblyRear.torque.value

        dv = cogForce ./ mass.value - angVel × vel
        dangularVelocity = cogMoment ./ inertia.value
        dx = [dv[1], dv[2], angVel[3], dangularVelocity[3]]
        return dx
    end

    state_descriptor = VarEntry[
        VarEntry("vx",    [velocity => 1],          :state),
        VarEntry("vy",    [velocity => 2], -1.5, 1.5, :state),
        VarEntry("psi",   [psi => 0],               :state),
        VarEntry("omega", [angularVelocity => 3],   :state),
        VarEntry("n",     [n => 0],                 :state),
        VarEntry("t",     [s => 0],                 :state),
    ]

    control_descriptor = VarEntry[
        VarEntry("torque_front", [drivetrain.motors[1].torque => 0], :control),
        #VarEntry("torque_rear",  [drivetrain.motors[2].torque => 0], :control),
        VarEntry("steering",    [wheelAssemblyFront.steeringAngle => 0], :control),
        VarEntry("brake",       [brakeCommand => 0], :control),
    ]

    nControls.value = Float64(length(control_descriptor))
    nStates.value   = Float64(length(state_descriptor))

    function controlMapping(controls::AbstractVector)
        apply_mapping!(control_descriptor, controls)
        frontBrake = brakeCommand.value * brakeBias.value
        rearBrake = brakeCommand.value * (1 - brakeBias.value)
        drivetrain.tires[1].brakingForce.value = frontBrake
        drivetrain.tires[2].brakingForce.value = rearBrake
    end

    function stateMapping(states::AbstractVector)
        apply_mapping!(state_descriptor, states)
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
        simplestSingleTrack,
        p,
        controlMapping,
        stateMapping,
        f -> (0.0),
        drivetrain,
        aero,
        suspension,
        chassis,
        [wheelAssemblyFront, wheelAssemblyRear],
        state_descriptor,
        control_descriptor,
    )
    return afto
end
