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

    function calculate()
        
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

    function calculate()
        totalWeight = chassis.mass.value
        frontRatio = 1 - chassis.CoG_X_pos.value
        rearRatio = chassis.CoG_X_pos.value
        leftRatio = 1 - chassis.CoG_Y_pos.value
        rightRatio = chassis.CoG_Y_pos.value
        Fz_FL = frontRatio * leftRatio * totalWeight
        Fz_FR = frontRatio * rightRatio * totalWeight
        Fz_RL = rearRatio * leftRatio * totalWeight
        Fz_RR = rearRatio * rightRatio * totalWeight
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