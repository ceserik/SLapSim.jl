using Revise
using SLapSim





function createTwintrack()

    gearboxFL = createCTU25gearbox()
    gearboxFR = createCTU25gearbox()
    gearboxRL = createCTU25gearbox()
    gearboxRR = createCTU25gearbox()

    motorFL = createFischerMotor()
    motorFR = createFischerMotor()
    motorRL = createFischerMotor()
    motorRR = createFischerMotor()

    tireFL = createR20lin()
    tireFR = createR20lin()
    tireRL = createR20lin()
    tireRR = createR20lin()

    drivetrain = Drivetrain(
        [motorFL, motorFR, motorRL, motorRR],
        [gearboxFL, gearboxFR, gearboxRL, gearboxRR],
        [tireFL, tireFR, tireRL, tireRR],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createDummySuspension()

    chassis = createCTU25chassis()
    wheelAssemblyFL = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2, chassis.track.value / 2, 0]))
    wheelAssemblyFR = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2, -chassis.track.value / 2, 0]))
    wheelAssemblyRL = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2, chassis.track.value / 2, 0]))
    wheelAssemblyRR = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2, -chassis.track.value / 2, 0]))

    wheelAssemblyFL.tire = tireFL
    wheelAssemblyFL.motor = motorFL
    wheelAssemblyFL.gearbox = gearboxFL
    wheelAssemblyFR.tire = tireFR
    wheelAssemblyFR.motor = motorFR
    wheelAssemblyFR.gearbox = gearboxFR
    wheelAssemblyRL.tire = tireRL
    wheelAssemblyRL.motor = motorRL
    wheelAssemblyRL.gearbox = gearboxRL
    wheelAssemblyRR.tire = tireRR
    wheelAssemblyRR.motor = motorRR
    wheelAssemblyRR.gearbox = gearboxRR

    for wa in [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL, wheelAssemblyRR]
        wa.tire.maxForce.value = wa.motor.torqueSpeedFunction(0.0) * wa.gearbox.ratio.value / wa.tire.radius.value
    end

    velocity = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "translational velocity", "m/s")
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 1.0], "angular velocity", "rad/s")

    mass = carParameter{carVar}(280.0, "Mass", "kg")
    motorForce = carParameter{carVar}(1000.0, "motorForce", "N")
    lateralForce = carParameter{carVar}(0.0, "lateral Force", "N")
    CL = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    velocity = carParameter{Vector{carVar}}([15.0, 0.0, 0.0], "Speed X", "m/s")
    powerLimit = carParameter{carVar}(80000.0, "PowerLimit", "W")
    psi = carParameter{carVar}(0.0, "heading", "rad")
    n = carParameter{carVar}(0.0, "Distance from centerline", "m")
    nControls = carParameter{carVar}(3.0, "number of controlled parameters", "-")
    inertia = carParameter{carVar}(100.0, "Inertia", "kg*m^2")
    nStates = carParameter{carVar}(6.0, "number of car states", "-")
    s = carParameter{carVar}(1.0, "longitudinal position on track", "-")


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

    function anyTrack(track::Union{Track,Nothing}=nothing, k::Union{Int64,Nothing,Float64}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)
        steeringAngle = car.wheelAssemblies[1].steeringAngle.value
        velocity = car.carParameters.velocity.value
        angularVelocity = car.carParameters.angularVelocity.value
        tireFL = car.drivetrain.tires[1]
        tireFR = car.drivetrain.tires[2]



        # Transformation of velocities from cog to wheels
        for wa in car.wheelAssemblies
            wa.setPivotVelocity(angularVelocity, velocity)
        end

        #Steer the wheels
        for i in eachindex(car.wheelAssemblies)
            car.wheelAssemblies[i].setTireSpeeds(car.drivetrain.tires[i])
        end


        #gearing of forces from motor to tire
        for i in eachindex(car.drivetrain.gearboxes)
            car.drivetrain.gearboxes[i].torqueIn.value = car.drivetrain.motors[i].torque.value
            car.drivetrain.gearboxes[i].f()
        end

        # update tire maxForce from current motor speed
        for wa in car.wheelAssemblies
            if !isnothing(wa.motor) && !isnothing(wa.gearbox) && !isnothing(wa.tire)
                wa.tire.maxForce.value = wa.motor.torqueSpeedFunction(wa.motor.angularFrequency.value) * wa.gearbox.ratio.value / wa.tire.radius.value
            end
        end

        #create constraints for motor

        for i in eachindex(car.drivetrain.tires)
            car.drivetrain.tires[i].tireFunction(car.drivetrain.gearboxes[i].torqueOut.value, optiModel)
        end


        if !isnothing(optiModel)

            # TODO this has to be enfroced only for car parts which are controls
            # motor torque limit
            car.drivetrain.motors[3].constraints(torqueFR, optiModel)
            car.drivetrain.motors[4].constraints(torqueRear, optiModel)

            #steering angle
            car.wheelAssemblies[1].constraints(optiModel)
            #car.wheelAssemblies[2].constraints(optiModel)

            #hitbox
            car.chassis.hitbox(car.carParameters.n.value, track, optiModel)

            for tire in car.drivetrain.tires
                tire.tireConstraints(optiModel)
            end
        end
        #tire Function, currently same as in bachelors thesis
        # calculate Fz on tires


        car.drivetrain.tires[1].forces.value[3] = 0.25 * car.carParameters.mass.value * 9.81
        car.drivetrain.tires[2].forces.value[3] = 0.25 * car.carParameters.mass.value * 9.81
        car.drivetrain.tires[3].forces.value[3] = 0.25 * car.carParameters.mass.value * 9.81
        car.drivetrain.tires[4].forces.value[3] = 0.25 * car.carParameters.mass.value * 9.81



        for i in eachindex(car.wheelAssemblies)
            car.wheelAssemblies[i].setPivotForce(car.drivetrain.tires[i])
        end


        cogForce = zero(car.wheelAssemblies[1].forces.value)
        cogMoment = zero(car.wheelAssemblies[1].pivot2CoG(car.drivetrain.tires[1].forces.value))
        for i in eachindex(car.wheelAssemblies)
            cogForce = cogForce .+ car.wheelAssemblies[i].forces.value
            cogMoment = cogMoment + car.wheelAssemblies[i].pivot2CoG(car.drivetrain.tires[i].forces.value)
        end


        #enforce hitbox

        #print cogMoment and cogForce
        #println("CoG Moment 1: ", cogMoment1)
        #println("CoG Moment 2: ", cogMoment2)
        #println("CoG Force: ", cogForce)

        dv = cogForce / car.carParameters.mass.value + angularVelocity × velocity #really check what sign should be here !!!! podla mna bednarik skripta fyzika1 Kapitola 8  Neinerciální vztažné soustavy, neboli pro vyjádření časové změny libovolné vektorové veličiny v nečárkované soustavě je možné použít následujícího operátoru:  d·  dt = d′·  dt + ω × · , (8.11)
        dangularVelocity = cogMoment / car.carParameters.inertia.value
        dx = [dv[1], dv[2], angularVelocity[3], dangularVelocity[3]]
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
