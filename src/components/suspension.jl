mutable struct Suspension
    tlong::carParameter{carVar}
    tlat::carParameter{carVar}
    theave::carParameter{carVar}

    stiffnessLong::carParameter{carVar}
    dampingLong::carParameter{carVar}

    stiffnessLat::carParameter{carVar}
    dampingLat::carParameter{carVar}

    stiffnessHeave::carParameter{carVar}
    dampingHeave::carParameter{carVar}
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