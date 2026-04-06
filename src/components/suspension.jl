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
    lateralTransfer::carParameter{carVar}
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
    lateralTransfer = carParameter{carVar}(0.0,"lateral load transfer","N")

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
        lateralTransfer,
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
    lateralTransfer = carParameter{carVar}(0.0,"lateral load transfer","N")
    chassis::Union{Chassis, Nothing} = nothing


    function setInput(chassis_in::Chassis)
        chassis = chassis_in

    end

    function calculate(downforce::carVar=0.0, CoP::carVar=0.5)
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
        lateralTransfer,
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
    lateralTransfer = carParameter{carVar}(0.0,"lateral load transfer","N")
    chassis::Union{Chassis, Nothing} = nothing

    function setInput(chassis_in::Chassis)
        chassis = chassis_in
    end

    function calculate(downforce::carVar=0.0, CoP::carVar=0.5)
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
        lateralTransfer,
        calculate,
        setInput    )
    return susp
end


function createQuasi_steady_Suspension()
    tlong = carParameter{carVar}(0.0,"longitudinal transfer time constant","s")
    tlat = carParameter{carVar}(0.0,"lateral transfer time constant","s")
    theave = carParameter{carVar}(0.0,"heave transfer time constant","s")

    stiffnessLong = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    stiffnessLat = carParameter{carVar}(0.0,"Lat stifness","N/rad??")
    stiffnessHeave = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    
    dampingLong = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingLat = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingHeave = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    lateralTransfer = carParameter{carVar}(0.0,"lateral load transfer","N")
    chassis::Union{Chassis, Nothing} = nothing


    function setInput(chassis_in::Chassis)
        chassis = chassis_in

    end

    function calculate(ax,downforce,CoP)
        frontRatio = chassis.CoG_X_pos.value
        rearRatio = 1 - chassis.CoG_X_pos.value
        leftRatio = 1 - chassis.CoG_Y_pos.value
        rightRatio = chassis.CoG_Y_pos.value

        mass = chassis.mass.value
        weight = mass * 9.81
        weightFront = weight * frontRatio
        weightRear = weight * rearRatio

        aeroFront = downforce * CoP
        aeroRear = downforce * (1 - CoP)

        # longitudinal load transfer (N): m*ax*h/L, shifted from front to rear
        h_cog = chassis.CoG_Z_pos.value
        transfer = mass * ax * h_cog / chassis.wheelbase.value

        Fz_FL = (weightFront - transfer) * leftRatio  + aeroFront/2
        Fz_FR = (weightFront - transfer) * rightRatio + aeroFront/2
        Fz_RL = (weightRear  + transfer) * leftRatio  + aeroRear/2
        Fz_RR = (weightRear  + transfer) * rightRatio + aeroRear/2

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
        lateralTransfer,
        calculate,
        setInput    )
    return susp
end