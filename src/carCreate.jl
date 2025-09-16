
using Interpolations

function lessContraint(a,b,model=nothing)
    if isnothing(model)
        a = min(a,b)
    else
        @constraint(model,a<=b)
    end
    return a
end

function greaterContraint(a,b,model=nothing)
    if isnothing(model)
        a = max(a,b)
    else
        @constraint(model,a>=b)
    end
    return a
end

function equalConstraint(a,b)
    if isnothing(model)
        a = b
    else
        @constraint(model,a==b)
    end    
    return a
end

function massPointCar(car,track,k, optiModel=nothing)
    # Simple  Mass point friction circle
    # inputs is a struct of all posible inputs of vehicle, controls are not separated from inputs
    # Get inputs
    m          = car.carParameters.mass.value
    inputForce = car.carParameters.motorForce.value
    CL         = car.carParameters.CL.value
    CD         = car.carParameters.CD.value
    maxPower   = car.carParameters.powerLimit.value

    # Get state
    vx         = car.carParameters.velocity.value[1]
        # this simplificqation is necceseary so that  universal function time -> path can be used
        #it means that car is always pointing the sma edirection as track and no sideways motion is possible
    car.carParameters.psi.value = track.theta[k]
    
    # Get track inputs
    c          = track.curvature[k]
    rho        = track.rho[k]
    μ          = track.μ[k]
    # Calculate forces

    Fz = 1/2 * rho * CL * vx^2
    Fy = m*vx^2*c

    maxMotorForce = 6507 #calculated for ctu25 should be added to inputs, vx torque characerisitic
    FxPowerMax = maxPower/(vx)
    #print(Fy)
    FxMaxsquared = max((Fz*μ + m*9.81*μ)^2 - Fy^2,0)
    # Add optimization constraints if model is provided

    if optiModel !== nothing
        ## tu urobit funkciu do ktorej dam obmedzenie a ona mi o spravi aby som nemusel stale davat ify
        ## aj ked to mozno je jedno lebo tento mass point bude mozno malo pouzivany?
        ## ale radsej to spravit, nech netreba prepisovat model lebo sa z toho zblaznim ked tam bude nieco inak
        @constraint(optiModel, (Fy/maxMotorForce)^2 + (inputForce/maxMotorForce)^2 <= ((Fz*μ + m*9.81*μ)/maxMotorForce)^2)
        @constraint(optiModel,inputForce/FxPowerMax <= FxPowerMax/FxPowerMax)
    else
        inputForce = min(inputForce, sqrt(FxMaxsquared))
        inputForce = max(inputForce, -sqrt(FxMaxsquared))
        inputForce = min(inputForce,maxMotorForce)
        inputForce = min(inputForce,FxPowerMax)
    end         

   # lessContraint((Fy/maxMotorForce)^2 + (inputForce/maxMotorForce)^2 ,((Fz*μ + m*9.81*μ)/maxMotorForce)^2,optiModel)
   # lessContraint(inputForce/FxPowerMax , FxPowerMax/FxPowerMax,optiModel)

        dstates =  [(inputForce -(1/2 *rho*CD*vx^2))/m 0  0    0     0 0]
    return dstates
end


#Simple mass point model
function createCTU25_1D()
    mass = carParameter(280.0, "Mass", "kg")
    inertia = carParameter(100.0, "Inertia", "kg*m^2")
    motorForce = carParameter(1000.0, "motorForce", "N")
    lateralForce = carParameter(0.0, "lateral Force", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    velocity = carParameter([15.0,0,0],"Speed X","m/s")
    angularVelocity = carParameter([0.0, 0.0, 0.0], "angular velocity", "rad/s");
    powerLimit = carParameter(80000.0,"PowerLimit","W")
    vy = carParameter(0.0,"Speed Y","m/s")
    psi = carParameter(0.0,"heading","rad")
    n = carParameter(0.0,"Distance from centerline","m")
    nControls = carParameter(1.0,"number of controlled parameters","-")

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

    function controlMapping(car,u)
        car.carParameters.motorForce.value = u[1]
        #paramcopy.CL.value = controls[2]
        
    end

    function mapping(car,u,x)
        car.carParameters.motorForce.value = u[1]
        car.carParameters.vx.value = x[1]
    end

    function stateMapping(car,x)
        car.carParameters.velocity.value = [x[1] 0 0]
    end
    car = Car(
        massPointCar,
        p,
        controlMapping,
        stateMapping,
        mapping
    )

    return car
end


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
