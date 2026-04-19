
mutable struct Gearbox{F,G}
    ratio::carParameter{carVar}
    torqueIn::carParameter{carVar}
    angularFrequencyIn::carParameter{carVar}
    torqueOut::carParameter{carVar}
    angularFrequencyOut::carParameter{carVar}
    efficiency::carParameter{carVar}
    compute::F
    setTorque::G
end

Base.show(io::IO, ::MIME"text/plain", obj::Gearbox) = prettyPrintComponent(io, obj)

function createCTU25gearbox(ratio_p::Float64 = 11.46)
    ratio = carParameter{carVar}(ratio_p,"Gearbox ratio","-")
    torqueIn = carParameter{carVar}(0.0,"torque on motor","-")
    angularFrequencyIn = carParameter{carVar}(0.0,"velocity on motor","rad/s")

    torqueOut = carParameter{carVar}(0.0,"torque on wheel","-")
    angularFrequencyOut = carParameter{carVar}(0.0,"velocity on wheel","rad/s")
    efficiency = carParameter{carVar}(1.0,"efficiency","-")

    function compute()
        torqueOut.value = torqueIn.value * ratio.value * efficiency.value
        angularFrequencyOut.value = angularFrequencyIn.value / ratio.value
    end

    function setTorque(torque::carVar)
        torqueIn.value = torque
    end
    
    gearbox = Gearbox(
        ratio,
        torqueIn,
        angularFrequencyIn,
        torqueOut,
        angularFrequencyOut,
        efficiency,
        compute,
        setTorque
    )
end
