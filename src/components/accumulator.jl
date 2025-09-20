mutable struct Accumulator
    capacity::carParameter{Float64}
    maxPower::carParameter{Float64}
    minPower::carParameter{Float64}
    voltage::carParameter{Float64}
    current::carParameter{Float64}
    SoC::carParameter{Float64}
    mass::carParameter{Float64}
    resistance::carParameter{Float64}
end
Base.show(io::IO, ::MIME"text/plain", obj::Accumulator) = prettyPrintComponent(io, obj)


function createPepikCTU25()
    capacity = carParameter(7.16,"ACP capacity","kWh")
    maxPower = carParameter(80.0,"max Power out","kW")
    minPower = carParameter(80.0,"min Power out","kW")
    voltage  = carParameter(600.0,"Voltage","kW")
    current  = carParameter(0.0,"Current","kW")
    SoC = carParameter(100.0,"State of Charge","%")
    mass = carParameter(38.49,"mass","kg")
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
