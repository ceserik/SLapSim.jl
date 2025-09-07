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

Base.show(io::IO, ::MIME"text/plain", obj::Suspension) = prettyPrintComponent(io, obj)

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