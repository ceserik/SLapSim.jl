using Revise
using SLapSim


function createSimplestSingleTrack()

    gearboxFront = createCTU25gearbox()
    gearboxRear = createCTU25gearbox()

    motorFront = createFischerMotor()
    motorRear = createFischerMotor()
    
    tireFront = createR20lin(motorFront,gearboxFront)
    tireRear = createR20lin(motorRear,gearboxRear)


    drivetrain = Drivetrain(
        [motorFront, motorRear],
        [gearboxFront, gearboxRear],
        [tireFront, tireRear],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createDummySuspension()
    chassis = createCTU25chassis()
    wheelAssemblyFront = createBasicWheelAssembly(Vector{carVar}([1.520 / 2, 0, 0]))
    wheelAssemblyRear = createBasicWheelAssembly(Vector{carVar}([-1.520 / 2, 0, 0]))


    velocity = carParameter{Vector{carVar}}([15.0, 0.0, 0.0], "Speed X", "m/s")
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 1.0], "angular velocity", "rad/s")

    mass = carParameter{carVar}(280.0, "Mass", "kg")
    motorForce = carParameter{carVar}(1000.0, "motorForce", "N")
    lateralForce = carParameter{carVar}(0.0, "lateral Force", "N")
    CL = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    powerLimit = carParameter{carVar}(80000.0, "PowerLimit", "W")
    psi = carParameter{carVar}(0.0, "heading", "rad")
    n = carParameter{carVar}(0.0, "Distance from centerline", "m")
    nControls = carParameter{carVar}(3.0, "number of controlled parameters", "-")
    inertia = carParameter{carVar}(100.0, "Inertia", "kg*m^2")
    nStates = carParameter{carVar}(6.0, "number of car states", "-")
    s = carParameter{carVar}(1.0, "longitudinal position on track", "-")


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

        # calculate Fz on tires
        drivetrain.tires[1].forces.value[3] = 0.5 * mass.value * 9.81
        drivetrain.tires[2].forces.value[3] = 0.5 * mass.value * 9.81
        
        #tire forces
        drivetrain.tires[1].compute(gbFront.torqueOut.value, optiModel)
        drivetrain.tires[2].compute(gbRear.torqueOut.value, optiModel)

        if !isnothing(optiModel)
            drivetrain.motors[1].constraints(drivetrain.motors[1].torque.value, optiModel)
            drivetrain.motors[2].constraints(drivetrain.motors[2].torque.value, optiModel)
            wheelAssemblyFront.constraints(optiModel)
            wheelAssemblyRear.constraints(optiModel)
            chassis.hitbox(n.value, track, optiModel)
            drivetrain.tires[1].tireConstraints(optiModel)
            drivetrain.tires[2].tireConstraints(optiModel)
        end



        #propagate forces from tires to wheel assemblies to COG
#        @infiltrate
        wheelAssemblyFront.getTorque(drivetrain.tires[1].forces.value)
        wheelAssemblyRear.getTorque(drivetrain.tires[2].forces.value)

        cogForce = wheelAssemblyFront.forces.value .+ wheelAssemblyRear.forces.value
        #@infiltrate
        cogMoment = wheelAssemblyFront.torque.value + wheelAssemblyRear.torque.value

        dv = cogForce / mass.value - angVel × vel
        dangularVelocity = cogMoment / inertia.value
        dx = [dv[1], dv[2], angVel[3], dangularVelocity[3]]
        return dx
    end


    function controlMapping(controls)
        drivetrain.motors[1].torque.value = controls[1]
        drivetrain.motors[2].torque.value = controls[2]
        wheelAssemblyFront.steeringAngle.value = controls[3]
    end

    function stateMapping(states)
        velocity.value = [states[1], states[2], 0.0]
        psi.value = states[3]
        angularVelocity.value = [0.0, 0.0, states[4]]
        n.value = states[5]
        s.value = states[6]
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
    )
    return afto
end
