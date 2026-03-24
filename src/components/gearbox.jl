
mutable struct Gearbox{F,G}
    ratio::carParameter{carVar}
    torqueIn::carParameter{carVar}
    angularFrequencyIn::carParameter{carVar}
    torqueOut::carParameter{carVar}
    angularFrequencyOut::carParameter{carVar}
    loss::carParameter{carVar}
    compute::F
    setTorque::G
end

Base.show(io::IO, ::MIME"text/plain", obj::Gearbox) = prettyPrintComponent(io, obj)

function createCTU25gearbox()
    ratio = carParameter{carVar}(11.46,"Gearbox ratio","-")
    torqueIn = carParameter{carVar}(0.0,"torque in","-")
    speedIn = carParameter{carVar}(0.0,"velocity in","rad/s")

    torqueOut = carParameter{carVar}(0.0,"torque out","-")
    angularFrequencyOut = carParameter{carVar}(0.0,"velocity out","rad/s")
    loss = carParameter{carVar}(0.0,"loss","Watt")

    function compute()
        torqueOut.value = torqueIn.value * ratio.value
        angularFrequencyOut.value = speedIn.value / ratio.value
    end

    function setTorque(torque::carVar)
        torqueIn.value = torque
    end
    
    gearbox = Gearbox(
        ratio,
        torqueIn,
        speedIn,
        torqueOut,
        angularFrequencyOut,
        loss,
        compute,
        setTorque
    )
end
