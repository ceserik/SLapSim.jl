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



function createPepikCTU25()
    capacity = carParameter(7.16,"ACP capacity","kWh")
    maxPower = carParameter(80,"max Power out","kW")
    minPower = carParameter(80,"min Power out","kW")
    voltage  = carParameter(600,"Voltage","kW")
    current  = carParameter(0,"Current","kW")
    SoC = carParameter(100,"State of Charge","%")
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
