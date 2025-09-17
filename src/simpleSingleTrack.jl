function simplestSingleTrack(car,track=nothing,k=nothing, optiModel=nothing)
    # assign names for easier reading
    torqueFront = car.drivetrain.motors[1].torque.value
    torqueRear = car.drivetrain.motors[2].torque.value
    steeringAngle = car.wheelAssemblies[1].steeringAngle.value

    gbFront = car.drivetrain.gearboxes[1]
    gbRear = car.drivetrain.gearboxes[2]

    velocity = car.carParameters.velocity.value
    angularVelocity = car.carParameters.angularVelocity.value

    tireFront = car.drivetrain.tires[1]
    tireRear = car.drivetrain.tires[2]

    # Transformation of velocities from cog to wheels
    car.wheelAssemblies[1].setPivotVelocity(angularVelocity,velocity)
    car.wheelAssemblies[2].setPivotVelocity(angularVelocity,velocity)

    #Steer the wheels
    car.wheelAssemblies[1].setTireSpeeds(tireFront)
    car.wheelAssemblies[2].setTireSpeeds(tireRear)

    #gearing of forces from motor to tire, would be nice o have this in a loop
    gbFront.torqueIn.value = torqueFront
    gbFront.f()

    gbRear.torqueIn.value = torqueRear
    gbRear.f()

    #create constraints for motor
    if !isnothing(optiModel)
        car.drivetrain.motors[1].constrains(optiModel)
        car.drivetrain.motors[2].constrains(optiModel)
    end
    #tire Function, currently same as in bachelors thesis
    # calculate Fz on tires
    car.drivetrain.tires[1].forces.value[3] = 0.5 * car.carParameters.mass.value * 9.81
    car.drivetrain.tires[2].forces.value[3] = 0.5 * car.carParameters.mass.value * 9.81


    tireFront.tireFunction(gbFront.torqueOut.value)
    tireRear.tireFunction(gbRear.torqueOut.value)

    car.wheelAssemblies[1].setPivotForce(tireFront)
    car.wheelAssemblies[2].setPivotForce(tireRear)


    cogMoment1 = car.wheelAssemblies[1].pivot2CoG(tireFront.forces.value) # zero at the end is because vertcial force would cause car to spn around y
    cogMoment2 = car.wheelAssemblies[2].pivot2CoG(tireRear.forces.value)

    cogForce = car.wheelAssemblies[1].forces.value .+ car.wheelAssemblies[2].forces.value
    cogMoment = cogMoment1 + cogMoment2

    #print cogMoment and cogForce
    #println("CoG Moment 1: ", cogMoment1)
    #println("CoG Moment 2: ", cogMoment2)
    #println("CoG Force: ", cogForce)

    dv = cogForce/car.carParameters.mass.value + angularVelocity × velocity
    dangularVelocity = cogMoment/car.carParameters.inertia.value
    dx = [dv[1],dv[2],angularVelocity[3],dangularVelocity[3]]
    return dx

end


function createSimplestSingleTrack()

    gearboxFront = createCTU25gearbox()
    gearboxRear = createCTU25gearbox()

    motorFront = createFischerMotor()
    motorRear = createFischerMotor()

    #get max motor torque for scaling
    maxMotorTorqueFront = motorFront.torqueSpeedFunction(0) * gearboxFront.ratio.value
    maxMotorTorqueRear = motorFront.torqueSpeedFunction(0) * gearboxRear.ratio.value

    tireFront = createR20lin(maxMotorTorqueFront)
    tireRear = createR20lin(maxMotorTorqueRear)


    drivetrain = Drivetrain(
        [motorFront,motorRear],
        [gearboxFront, gearboxRear],
        [tireFront,tireRear],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createDummySuspension()
    wheelAssemblyFront = createBasicWheelAssembly([1.520/2, 0, 0]) # wheelbase musi byt parameter!!!!
    wheelAssemblyRear = createBasicWheelAssembly([-1.520/2, 0, 0])
    chassis = createCTU25chassis()

    velocity = carParameter([10.0, 10.0, 0.0], "translational velocity", "m/s");
    angularVelocity = carParameter([0.0, 0.0, 1.0], "angular velocity", "rad/s");
    
    mass = carParameter(280.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    lateralForce = carParameter(0.0, "lateral Force", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    velocity = carParameter([15,0,0],"Speed X","m/s")
    powerLimit = carParameter(80000.0,"PowerLimit","W")
    psi = carParameter(0.0,"heading","rad")
    n = carParameter(0.0,"Distance from centerline","m")
    nControls = carParameter(2.0,"number of controlled parameters","-")
    inertia = carParameter(100.0, "Inertia", "kg*m^2")


    function controlMapping(car, controls)
        car.drivetrain.motors[1].torque.value = controls[1]  
        car.drivetrain.motors[2].torque.value = controls[1]  
        return car
    end
    function stateMapping(car,states)
        car.carParameters.velocity.value = [states[1], states[2], 0.0]
        car.carParameters.psi.value = states[3]
        car.carParameters.angularVelocity.value = [0.0, 0.0, states[4]]
        return car
    end
    
    p = carParameters(
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
        nControls
    )

    afto = Car2(
        drivetrain,
        aero,
        suspension,
        chassis,
        [wheelAssemblyFront,wheelAssemblyRear],
        p,
        simplestSingleTrack,
        controlMapping,
        stateMapping
    )
    return afto
end
