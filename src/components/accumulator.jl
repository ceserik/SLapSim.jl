mutable struct Accumulator{F1}
    capacity::carParameter{carVar}
    power::carParameter{carVar}
    voltage::carParameter{carVar}
    SoC::carParameter{carVar}
    mass::carParameter{carVar}
    resistance::carParameter{carVar}
    compute::F1
end
Base.show(io::IO, ::MIME"text/plain", obj::Accumulator) = prettyPrintComponent(io, obj)


function createAccumulator()
    capacity   = carParameter{carVar}(7.16, "ACP capacity", "kWh")
    power      = carParameter{carVar}(0.0, "Power", "W", :static, [-80_000.0, 80_000.0])
    voltage    = carParameter{carVar}(600.0, "Voltage", "V")
    SoC        = carParameter{carVar}(100.0, "State of Charge", "%")
    mass       = carParameter{carVar}(38.49, "mass", "kg")
    resistance = carParameter{carVar}(0.01, "Internal Resistance", "ohm")

    function compute(motors, model = nothing)
        # motor.computeElectric must be called beforehand to populate powerElectrical
        power.value = sum(m.powerElectrical.value for m in motors)
        if !isnothing(model)
            sf = power.limits[2]
            @constraint(model, power.value / sf <= power.limits[2] / sf)
            @constraint(model, power.value / sf >= power.limits[1] / sf)
        end
    end

    Accumulator(
        capacity,
        power,
        voltage,
        SoC,
        mass,
        resistance,
        compute
    )
end
