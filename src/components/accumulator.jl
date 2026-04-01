mutable struct Accumulator
    capacity::carParameter{carVar}
    maxPower::carParameter{carVar}
    minPower::carParameter{carVar}
    voltage::carParameter{carVar}
    current::carParameter{carVar}
    SoC::carParameter{carVar}
    mass::carParameter{carVar}
    resistance::carParameter{carVar}
end
Base.show(io::IO, ::MIME"text/plain", obj::Accumulator) = prettyPrintComponent(io, obj)


function createPepikCTU25()
    capacity = carParameter{carVar}(7.16,"ACP capacity","kWh")
    maxPower = carParameter{carVar}(80.0,"max Power out","kW")
    minPower = carParameter{carVar}(80.0,"min Power out","kW")
    voltage  = carParameter{carVar}(600.0,"Voltage","kW")
    current  = carParameter{carVar}(0.0,"Current","kW")
    SoC = carParameter{carVar}(100.0,"State of Charge","%")
    mass = carParameter{carVar}(38.49,"mass","kg")
    resistance = carParameter{carVar}(0.01,"Internal Resistance???","ohm")

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
