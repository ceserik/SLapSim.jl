# Define a Parameter struct
using LinearAlgebra

mutable struct carParameter
    value::Any
    name::String
    unit::String
    # Constructor for numeric values
    #carParameter(value::T, name::String, unit::String) where T = new(value, name, unit)
end

mutable struct carParameters
    mass::carParameter
    motorForce::carParameter
    CL::carParameter
    CD::carParameter
    vx::carParameter
    vy::carParameter
    psi::carParameter
    n::carParameter
    powerLimit::carParameter
    lateralForce::carParameter
    nControls::carParameter
end

mutable struct Tire
    radius::carParameter
    width::carParameter
    inertia::carParameter
    mass::carParameter
    angularVelocity::carParameter
    vertForce::carParameter
    longForce::carParameter
    latForce::carParameter
    slipAngle::carParameter
    slipRatio::carParameter
    tireFunction
end

# Constructor for Tire
function Tire(
    radius::carParameter,
    width::carParameter,
    inertia::carParameter,
    mass::carParameter,
    angularVelocity::carParameter,
    vertForce::carParameter,
    longForce::carParameter,
    latForce::carParameter,
    slipAngle::carParameter,
    slipRatio::carParameter,
    tireFunction)
    # where T
    #ire(radius, width, inertia, mass, angularVelocity, vertForce, longForce, latForce, slipAngle, slipRatio, tireFunction)
end

mutable struct Drivetrain
    motors
    gearboxes
    tires
    accumulators
end

mutable struct Motor
    torque::carParameter
    angularVelocity::carParameter
    mass::carParameter
    loss::carParameter
    torqueSpeedFunction #mapping speed to max torque
end

mutable struct Gearbox
    ratio::carParameter
    torqueIn::carParameter
    angularVelocityIn::carParameter
    torqueOut::carParameter
    angularVelocityOut::carParameter
    loss::carParameter
    gearboxFunction
    
end

mutable struct Accumulator
    capacity::carParameter
    maxPower::carParameter
    minPower::carParameter
    voltage::carParameter
    current::carParameter
    SoC::carParameter
    mass::carParameter
    resistance::carParameter
end

mutable struct Aero
    CL::carParameter
    CD::carParameter
    CoP::carParameter
end

mutable struct Suspension
    tlong::carParameter
    tlat::carParameter
    theave::carParameter

    stiffnessLong::carParameter
    dampingLong::carParameter

    stiffnessLat::carParameter
    dampingLat::carParameter

    stiffnessHeave::carParameter
    dampingHeave::carParameter
    # da sa posuvat pitch center???
end


mutable struct WheelAssembly
    position
    velocity
    CoG2wheelAssembly
    wheel2CoG
    steeringAngle::carParameter
    rotZ
end



# Define the Car struct with parameters
mutable struct Car
    carFunction
    carParameters
    controlMapping
    stateMapping
    mapping

    # Inner constructor
    function Car(carFunction, carParameters, controlMapping=nothing, stateMapping=nothing, mapping=nothing)
        new(carFunction, carParameters, controlMapping, stateMapping, mapping)
    end
end



mutable struct Car2
    drivetrain
    aero
    suspension
    wheelAssemblies

end



function createDummySuspension()
    tlong = carParameter(0.0,"longitudinal transfer time constant","s")
    tlat = carParameter(0.0,"lateral transfer time constant","s")
    theave = carParameter(0.0,"heave transfer time constant","s")

    stiffnessLong = carParameter(0.0,"Long stifness","N/rad??")
    stiffnessLat = carParameter(0.0,"Lat stifness","N/rad??")
    stiffnessHeave = carParameter(0.0,"Long stifness","N/rad??")
    
    dampingLong = carParameter(0.0,"long dampnng","N/rad/s`")
    dampingLat = carParameter(0.0,"long dampnng","N/rad/s`")
    dampingHeave = carParameter(0.0,"long dampnng","N/rad/s`")

    susp = Suspension(
        tlong,
        tlat,
        theave,
        stiffnessLong,
        dampingLong,
        stiffnessLat,
        dampingLat,
        stiffnessHeave,
        dampingHeave    )
    return susp
end

function createBasicAero()
    CL = carParameter(-5.0,"Lift coeffcient","-")
    CD = carParameter(2,"Drag coeffcient","-")
    CoP = carParameter(0.5,"Cener of pressure on front","-")
    aero = Aero(
        CL,
        CD,
        CoP
    )
    return aero

end

function createBasicWheelAssembly(position)
    steeringAngle = carParameter(0.0,"steering angle","rad")

    function rotZ(angle)
        out = [
            cos(angle) -sin(angle) 0;
            sin(angle)  cos(angle) 0;
            0           0          1
        ]
        return out
    end

    function CoG2wheelAssembly(position, velocity, angularVelocity)
        #returns speed at pivot point
        return velocity + cross(angularVelocity, position)
    end
    function wheel2CoG(forces, position)
        #returns moments on cog from wheel forces
        moments = cross(position, forces)
        return moments
    end

    testWheelAssembly = WheelAssembly(
        position,
        [0 0 0],
        CoG2wheelAssembly,
        wheel2CoG,
        steeringAngle,
        rotZ
    )
    return testWheelAssembly
end


function createPepikCTU25()
    capacity = carParameter(7.1,"ACP capacity","kWh")
    maxPower = carParameter(80,"max Power out","kW")
    minPower = carParameter(80,"min Power out","kW")
    voltage  = carParameter(600,"Voltage","kW")
    current  = carParameter(0,"Current","kW")
    SoC = carParameter(100,"State of Charge","%")
    mass = carParameter(44,"mass???","kg")
    resistance = carParameter(0.01,"Internal Resistance???","ohm")

    ACP = Accumulator(
        capacity,
        maxPower,
        minPower,
        voltage,
        current,
        SoC,
        mass,
        resistance
    )
    return ACP
end

function createR20lin()
    radius = carParameter(0.205, "tire radius", "m")
    width = carParameter(0.3, "tire width, wrong", "m")
    inertia = carParameter(0.3, "tire width, wrong", "m")
    mass = carParameter(1.0, "tire mass,wrong", "kg")
    angularVelocity = carParameter(0.0, "angular velocity", "rad/s")
    vertForce = carParameter(0.0, "Vertical force on tire", "N")
    longForce = carParameter(0.0, "Vertical force from tire", "N")
    latForce = carParameter(0.0, "Vertical force from tire", "N")
    slipAngle = carParameter(0.0, "slip angle", "rad")
    slipRatio = carParameter(0.0, "slip ratio", "-")

    function tireFunction(tire, velocity, vertForce)
        tire.slipAngle = atan(velocity[1], velocity[2])
        tire.slipRatio = tire.angularVelocity * tire.radius / velocity[1]
        tire.latForce = tire.slipAngle * vertForce
        tire.longForce = tire.slipRatio * vertForce
    end

    return Tire(
        radius,
        width,
        inertia,
        mass,
        angularVelocity,
        vertForce,
        longForce,
        latForce,
        slipAngle,
        slipRatio,
        tireFunction
    )
end

function createFischerMotor()
    torque = carParameter(0.0,"motor torque","Nm")
    angularVelocity = carParameter(0.0,"angular velocity","rad/s")
    loss = carParameter(0.0,"loss","W")
    torqueSpeedFunction = f(angularVelocity) = 29
    mass = carParameter(3.0,"mass??","kg")
    motor = Motor(
        torque,
        angularVelocity,
        mass,
        loss,
        torqueSpeedFunction
    )
    return motor
end

function createCTU25gearobx()
    ratio = carParameter(11.46,"Gearbox ratio","-")
    torqueIn = carParameter(0,"torque in","-")
    speedIn = carParameter(0,"velocity in","rad/s")

    torqueOut = carParameter(0,"torque out","-")
    speedOut = carParameter(0,"velocity out","rad/s")
    loss = carParameter(0,"loss","W")

    function gearboxFunction(gearbox) 
        gearbox.torqueOut = gearbox.torqueIn * 11.46
        gearbox.speedOut = gearbox.speedOut / 11/46
    end
    gearbox = Gearbox(
        ratio,
        torqueIn,
        speedIn,
        torqueOut,
        speedOut,
        loss,
        gearboxFunction
    )
end

function createCTU25_2D()
    mass = carParameter(280.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    lateralForce = carParameter(0.0, "lateral Force", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    vx = carParameter(15.0,"Speed X","m/s")
    powerLimit = carParameter(80000.0,"PowerLimit","W")
    vy = carParameter(0.0,"Speed Y","m/s")
    psi = carParameter(0.0,"heading","rad")
    n = carParameter(0.0,"Distance from centerline","m")
    nControls = carParameter(2.0,"number of controlled parameters","-")
    p = carParameters(
        mass,
        motorForce,
        CL,
        CD,
        vx,
        vy,
        psi,
        n,
        powerLimit,
        lateralForce,
        nControls
    )

    function controlMapping(input, controls)
        input.motorForce.value = controls[1]  
        #paramcopy.CL.value = controls[2]
        return input
    end

    function mapping(car,u,x)
        car.carParameters.motorForce.value = u[1]
        car.carParameters.lateralForce.value = u[2]

        car.carParameters.vx.value = x[1]
        car.carParameters.vy.value = x[2]

    end

    function stateMapping(input,states)
        input.vx.value = states[1]
        #map psi to track heading here, maybe input should be also track mapping should be unified
        return input
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
    


function createSimplestSingleTrack()
    tireFront = createR20lin()
    tireRear = createR20lin()
    gearboxFront = createCTU25gearobx()
    gearboxRear = createCTU25gearobx()

    motorFront = createFischerMotor()
    motorRear = createFischerMotor()

    drivetrain = Drivetrain(
        [motorFront,motorRear],
        [gearboxFront, gearboxRear],
        [tireFront,tireRear],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createDummySuspension()
    wheelAssemblyFront = createBasicWheelAssembly([1.520/2 0 0])
    wheelAssemblyRear = createBasicWheelAssembly([-1.520/2 0 0])
    afto = Car2(
        drivetrain,
        aero,
        suspension,
        [wheelAssemblyFront,wheelAssemblyRear]
    )

end



print("carParams.jl great success")
