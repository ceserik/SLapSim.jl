using Revise
using SLapSim


function simpleTwinTrack(
    car::Car,
    track::Union{Track,Nothing}=nothing,
    k::Union{Int64,Nothing,Float64}=nothing,
    optiModel::Union{JuMP.Model,Nothing}=nothing)
    # assign names for easier reading
    T_FL = car.drivetrain.motors[1].torque.value
    T_FR = car.drivetrain.motors[2].torque.value
    T_RL = car.drivetrain.motors[3].torque.value
    T_RR = car.drivetrain.motors[4].torque.value
    steeringAngle = car.wheelAssemblies[1].steeringAngle.value

    gbFL = car.drivetrain.gearboxes[1]
    gbFR = car.drivetrain.gearboxes[2]
    gbRL = car.drivetrain.gearboxes[3]
    gbRR = car.drivetrain.gearboxes[4]

    velocity = car.carParameters.velocity.value
    angularVelocity = car.carParameters.angularVelocity.value

    tireFL = car.drivetrain.tires[1]
    tireFR = car.drivetrain.tires[2]
    tireRL = car.drivetrain.tires[3]
    tireRR = car.drivetrain.tires[4]



    # Transformation of velocities from cog to wheels
    car.wheelAssemblies[1].setPivotVelocity(angularVelocity,velocity)
    car.wheelAssemblies[2].setPivotVelocity(angularVelocity,velocity)
    car.wheelAssemblies[3].setPivotVelocity(angularVelocity,velocity)
    car.wheelAssemblies[4].setPivotVelocity(angularVelocity,velocity)

    #Steer the wheels
    car.wheelAssemblies[1].setTireSpeeds(tireFL)
    car.wheelAssemblies[2].setTireSpeeds(tireFR)
    car.wheelAssemblies[3].setTireSpeeds(tireRL)
    car.wheelAssemblies[4].setTireSpeeds(tireRR)


    #gearing of forces from motor to tire, would be nice o have this in a loop
    
    gbFL.torqueIn.value = T_FL
    gbFL.f()

    gbFR.torqueIn.value = T_FR
    gbFR.f()

    gbRL.torqueIn.value = T_RL
    gbRL.f()

    gbFL.torqueIn.value = T_RR
    gbRR.f()

    #create constraints for motor
    tireFL.tireFunction(gbFL.torqueOut.value,optiModel)
    tireFR.tireFunction(gbFR.torqueOut.value,optiModel)
    tireRL.tireFunction(gbRL.torqueOut.value,optiModel)
    tireRR.tireFunction(gbRR.torqueOut.value,optiModel)

    
    if !isnothing(optiModel)
        # motor torque limit
        car.drivetrain.motors[3].constraints(torqueFR,optiModel)
        car.drivetrain.motors[4].constraints(torqueRear,optiModel)

        #steering angle
        car.wheelAssemblies[1].constraints(optiModel)
        car.wheelAssemblies[2].constraints(optiModel)
        

        #hitbox
        car.chassis.hitbox(car.carParameters.n.value,track,optiModel)


        tireFront.tireConstraints(optiModel)
        tireRear.tireConstraints(optiModel)
    end
    #tire Function, currently same as in bachelors thesis
    # calculate Fz on tires
    

    car.drivetrain.tires[1].forces.value[3] = 0.5 * car.carParameters.mass.value * 9.81
    car.drivetrain.tires[2].forces.value[3] = 0.5 * car.carParameters.mass.value * 9.81


    

    car.wheelAssemblies[1].setPivotForce(tireFront)
    car.wheelAssemblies[2].setPivotForce(tireRear)


    cogMoment1 = car.wheelAssemblies[1].pivot2CoG(tireFront.forces.value) # zero at the end is because vertcial force would cause car to spn around y
    cogMoment2 = car.wheelAssemblies[2].pivot2CoG(tireRear.forces.value)

    cogForce = car.wheelAssemblies[1].forces.value .+ car.wheelAssemblies[2].forces.value
    cogMoment = cogMoment1 + cogMoment2


    #enforce hitbox
    
    #print cogMoment and cogForce
    #println("CoG Moment 1: ", cogMoment1)
    #println("CoG Moment 2: ", cogMoment2)
    #println("CoG Force: ", cogForce)

    dv = cogForce/car.carParameters.mass.value + angularVelocity × velocity #really check what sign should be here !!!! podla mna bednarik skripta fyzika1 Kapitola 8  Neinerciální vztažné soustavy, neboli pro vyjádření časové změny libovolné vektorové veličiny v nečárkované soustavě je možné použít následujícího operátoru:  d·  dt = d′·  dt + ω × · , (8.11)
    dangularVelocity = cogMoment/car.carParameters.inertia.value
    dx = [dv[1],dv[2],angularVelocity[3],dangularVelocity[3]]
    return dx

end


function createSimplestSingleTrack()

    gearboxFL = createCTU25gearbox()
    gearboxFR = createCTU25gearbox()
    gearboxRL = createCTU25gearbox()
    gearboxRR = createCTU25gearbox()

    motorFL = createFischerMotor()
    motorFR = createFischerMotor()
    motorRL = createFischerMotor()
    motorRR = createFischerMotor()

    #get max motor torque for scaling
    maxMotorTorqueFL = motorFL.torqueSpeedFunction(0.0) * gearboxFL.ratio.value
    maxMotorTorqueFR = motorFR.torqueSpeedFunction(0.0) * gearboxFR.ratio.value

    maxMotorTorqueRL = motorRL.torqueSpeedFunction(0.0) * gearboxRL.ratio.value
    maxMotorTorqueRR = motorRR.torqueSpeedFunction(0.0) * gearboxRR.ratio.value
    

    
    tireFL = createR20lin(maxMotorTorqueFL)
    tireFR = createR20lin(maxMotorTorqueFR)
    tireRL = createR20lin(maxMotorTorqueRL)
    tireRR = createR20lin(maxMotorTorqueRR)



    drivetrain = Drivetrain(
        [motorFL,motorFR,motorRL,motorRR],
        [gearboxFL, gearboxFR,gearboxRL,gearboxRR],
        [tireFL,tireFR,tireRL,tireRR],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createDummySuspension()

    chassis = createCTU25chassis()
    wheelAssemblyFL = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value/2, chassis.track.value/2, 0])) # wheelbase musi byt parameter!!!!
    wheelAssemblyFR = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value/2, -chassis.track.value/2, 0])) # wheelbase musi byt parameter!!!!

    wheelAssemblyRL = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value/2, chassis.track.value/2, 0])) # wheelbase musi byt parameter!!!!
    wheelAssemblyRR = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value/2, -chassis.track.value/2, 0])) # wheelbase musi byt parameter!!!!


    
    

    velocity = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "translational velocity", "m/s");
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 1.0], "angular velocity", "rad/s");
    
    mass = carParameter{carVar}(280.0, "Mass", "kg")
    motorForce = carParameter{carVar}(1000.0, "motorForce", "N")
    lateralForce = carParameter{carVar}(0.0, "lateral Force", "N")
    CL = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    velocity = carParameter{Vector{carVar}}([15.0,0.0,0.0],"Speed X","m/s")
    powerLimit = carParameter{carVar}(80000.0,"PowerLimit","W")
    psi = carParameter{carVar}(0.0,"heading","rad")
    n = carParameter{carVar}(0.0,"Distance from centerline","m")
    nControls = carParameter{carVar}(3.0,"number of controlled parameters","-")
    inertia = carParameter{carVar}(100.0, "Inertia", "kg*m^2")
    nStates = carParameter{carVar}(6.0,"number of car states","-")
    s = carParameter{carVar}(1.0,"longitudinal position on track","-")


    function controlMapping(car, controls)
        car.drivetrain.motors[3].torque.value = controls[1]  
        car.drivetrain.motors[4].torque.value = controls[1]
        car.wheelAssemblies[1].steeringAngle.value = controls[2]
        car.wheelAssemblies[2].steeringAngle.value = controls[2]

    end
    function stateMapping(car, states)
        car.carParameters.velocity.value = [states[1], states[2], 0.0]
        car.carParameters.psi.value = states[3]
        car.carParameters.angularVelocity.value = [0.0, 0.0, states[4]]
        car.carParameters.n.value = states[5]
        car.carParameters.s.value = states[6]
    return car
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
        f->(0.0),
        drivetrain,
        aero,
        suspension,
        chassis,
        [wheelAssemblyFront,wheelAssemblyRear],
    )
    return afto
end
