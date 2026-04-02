mutable struct Suspension{F1,F2}
    tlong::carParameter{carVar}
    tlat::carParameter{carVar}
    theave::carParameter{carVar}

    stiffnessLong::carParameter{carVar}
    dampingLong::carParameter{carVar}

    stiffnessLat::carParameter{carVar}
    dampingLat::carParameter{carVar}

    stiffnessHeave::carParameter{carVar}
    dampingHeave::carParameter{carVar}
    calculate::F1
    setInput ::F2
    # da sa posuvat pitch center???
end

Base.show(io::IO, ::MIME"text/plain", obj::Suspension) = prettyPrintComponent(io, obj)

function createDummySuspension()
    tlong = carParameter{carVar}(0.0,"longitudinal transfer time constant","s")
    tlat = carParameter{carVar}(0.0,"lateral transfer time constant","s")
    theave = carParameter{carVar}(0.0,"heave transfer time constant","s")

    stiffnessLong = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    stiffnessLat = carParameter{carVar}(0.0,"Lat stifness","N/rad??")
    stiffnessHeave = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    
    dampingLong = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingLat = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingHeave = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")

    function setInput()
    
    
    end

    function calculate(downforce=0.0)

    end
    susp = Suspension(
        tlong,
        tlat,
        theave,
        stiffnessLong,
        dampingLong,
        stiffnessLat,
        dampingLat,
        stiffnessHeave,
        dampingHeave,
        calculate,
        setInput    )
    return susp
end




function createSimpleSuspension()
    tlong = carParameter{carVar}(0.0,"longitudinal transfer time constant","s")
    tlat = carParameter{carVar}(0.0,"lateral transfer time constant","s")
    theave = carParameter{carVar}(0.0,"heave transfer time constant","s")

    stiffnessLong = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    stiffnessLat = carParameter{carVar}(0.0,"Lat stifness","N/rad??")
    stiffnessHeave = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    
    dampingLong = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingLat = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingHeave = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    chassis::Union{Chassis, Nothing} = nothing


    function setInput(chassis_in::Chassis)
        chassis = chassis_in
    
    end

    function calculate(downforce=0.0, CoP=0.5)
        totalLoad = chassis.mass.value * 9.81 + downforce
        frontRatio = chassis.CoG_X_pos.value
        rearRatio = 1 - chassis.CoG_X_pos.value
        leftRatio = 1 - chassis.CoG_Y_pos.value
        rightRatio = chassis.CoG_Y_pos.value
        weightFront = chassis.mass.value * 9.81 * frontRatio
        weightRear = chassis.mass.value * 9.81 * rearRatio
        aeroFront = downforce * CoP
        aeroRear = downforce * (1 - CoP)
        Fz_FL = (weightFront + aeroFront) * leftRatio
        Fz_FR = (weightFront + aeroFront) * rightRatio
        Fz_RL = (weightRear + aeroRear) * leftRatio
        Fz_RR = (weightRear + aeroRear) * rightRatio
        return [Fz_FL, Fz_FR, Fz_RL, Fz_RR]
    end

    susp = Suspension(
        tlong,
        tlat,
        theave,
        stiffnessLong,
        dampingLong,
        stiffnessLat,
        dampingLat,
        stiffnessHeave,
        dampingHeave,
        calculate,
        setInput    )
    return susp
end

function createBusSuspension()
    tlong = carParameter{carVar}(0.0,"longitudinal transfer time constant","s")
    tlat = carParameter{carVar}(0.0,"lateral transfer time constant","s")
    theave = carParameter{carVar}(0.0,"heave transfer time constant","s")

    stiffnessLong = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    stiffnessLat = carParameter{carVar}(0.0,"Lat stifness","N/rad??")
    stiffnessHeave = carParameter{carVar}(0.0,"Long stifness","N/rad??")

    dampingLong = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingLat = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingHeave = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    chassis::Union{Chassis, Nothing} = nothing

    function setInput(chassis_in::Chassis)
        chassis = chassis_in
    end

    function calculate(downforce=0.0, CoP=0.5)
        frontRatio = chassis.CoG_X_pos.value
        rearRatio = 1 - chassis.CoG_X_pos.value
        leftRatio = 1 - chassis.CoG_Y_pos.value
        rightRatio = chassis.CoG_Y_pos.value
        weightFront = chassis.mass.value * 9.81 * frontRatio
        weightRear = chassis.mass.value * 9.81 * rearRatio
        aeroFront = downforce * CoP
        aeroRear = downforce * (1 - CoP)
        # front axle
        Fz_FL = (weightFront + aeroFront) * leftRatio
        Fz_FR = (weightFront + aeroFront) * rightRatio
        # rear: treat double axle as single, then split equally
        Fz_rear_left = (weightRear + aeroRear) * leftRatio / 2
        Fz_rear_right = (weightRear + aeroRear) * rightRatio / 2
        # order: FL, FR, RL1, RR1, RL2, RR2
        return [Fz_FL, Fz_FR, Fz_rear_left, Fz_rear_right, Fz_rear_left, Fz_rear_right]
    end

    susp = Suspension(
        tlong,
        tlat,
        theave,
        stiffnessLong,
        dampingLong,
        stiffnessLat,
        dampingLat,
        stiffnessHeave,
        dampingHeave,
        calculate,
        setInput    )
    return susp
end