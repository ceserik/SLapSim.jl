using Revise
using SLapSim
function createTwintrack()
    velocity = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s")
    angularVelocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s")

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

    gearboxFL = createCTU25gearbox()
    gearboxFR = createCTU25gearbox()
    gearboxRL = createCTU25gearbox()
    gearboxRR = createCTU25gearbox()

    motorFL = createFischerMotor()
    motorFR = createFischerMotor()
    motorRL = createFischerMotor()
    motorRR = createFischerMotor()

    tireFL = createR20lin(motorFL, gearboxFL) 
    tireFR = createR20lin(motorFR, gearboxFR)
    tireRL = createR20lin(motorRL, gearboxRL)
    tireRR = createR20lin(motorRR, gearboxRR)

    drivetrain = Drivetrain(
        [motorFL, motorFR, motorRL, motorRR],
        [gearboxFL, gearboxFR, gearboxRL, gearboxRR],
        [tireFL, tireFR, tireRL, tireRR],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createSimpleSuspension()

    chassis = createCTU25chassis()
    suspension.setInput(chassis)
    wheelAssemblyFL = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2, chassis.track.value / 2, 0]))
    wheelAssemblyFR = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2, -chassis.track.value / 2, 0]))
    wheelAssemblyRL = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2, chassis.track.value / 2, 0]))
    wheelAssemblyRR = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2, -chassis.track.value / 2, 0]))

    wheelAssemblies = [wheelAssemblyFL,wheelAssemblyFR,wheelAssemblyRL,wheelAssemblyRR]

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
        return car
    end

    function anyTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)

        #velocity = car.carParameters.velocity.value
        #angularVelocity = car.carParameters.angularVelocity.value
        # Transformation of velocities from cog to wheels and steering
        for wa in car.wheelAssemblies
            wa.setVelocity(angularVelocity.value, velocity.value)
        end

        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].setVelocity(wheelAssemblies[i].velocityTire.value)
        end

        #gearing of forces from motor to tire
        for i in eachindex(car.drivetrain.gearboxes)
            car.drivetrain.gearboxes[i].setTorque(car.drivetrain.motors[i].torque.value)
            car.drivetrain.gearboxes[i].compute()
        end



        rho = isnothing(track) ? RHO_SEA_LEVEL : track.rho[1]
        aeroForces = aero.compute(velocity.value[1], rho)

        forces = suspension.calculate(aeroForces.downforce, aero.CoP.value)
        for i in eachindex(forces)
            drivetrain.tires[i].forces.value[3] = forces[i]
        end

        #tire forces
        for i in eachindex(car.drivetrain.tires)
            car.drivetrain.tires[i].compute(car.drivetrain.gearboxes[i].torqueOut.value, optiModel)
        end

######################################################CONSTRAINTS###########################
        # TODO this has to be enfroced only for car parts which are controls
        # motor torque limit
        car.drivetrain.motors[3].constraints(drivetrain.motors[3].torque.value, optiModel)
        car.drivetrain.motors[4].constraints(drivetrain.motors[4].torque.value, optiModel)
        #steering angle
        car.wheelAssemblies[1].constraints(optiModel)
        #car.wheelAssemblies[2].constraints(optiModel)
        #hitbox
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

        dv = cogForce / car.carParameters.mass.value - angularVelocity.value × velocity.value #really check what sign should be here !!!! podla mna bednarik skripta fyzika1 Kapitola 8  Neinerciální vztažné soustavy, neboli pro vyjádření časové změny libovolné vektorové veličiny v nečárkované soustavě je možné použít následujícího operátoru:  d·  dt = d′·  dt + ω × · , (8.11)
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
