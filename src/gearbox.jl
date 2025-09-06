if !@isdefined(MODULES_INITIALIZED)
    include("initSlapSim.jl")
end

mutable struct Gearbox
    ratio::carParameter
    torqueIn::carParameter
    angularFrequencyIn::carParameter
    torqueOut::carParameter
    angularFrequencyOut::carParameter
    loss::carParameter
    f
end

# Pretty printing for Gearbox
import Base: show

function Base.show(io::IO, ::MIME"text/plain", obj::Gearbox)
    T = typeof(obj)
    println(io, "$(T)(")
    for field in fieldnames(T)
        println(io, "  $field = ", getfield(obj, field))
    end
    print(io, ")")
end


function createCTU25gearbox()
    ratio = carParameter(11.46,"Gearbox ratio","-")
    torqueIn = carParameter(0,"torque in","-")
    speedIn = carParameter(0,"velocity in","rad/s")

    torqueOut = carParameter(0,"torque out","-")
    angularFrequencyOut = carParameter(0,"velocity out","rad/s")
    loss = carParameter(0,"loss","Watt")

    function gearboxFunction() 
        gearbox.torqueOut.value = gearbox.torqueIn.value * 11.46
        gearbox.angularFrequencyOut.value = gearbox.angularFrequencyOut.value / 11/46
    end
    gearbox = Gearbox(
        ratio,
        torqueIn,
        speedIn,
        torqueOut,
        angularFrequencyOut,
        loss,
        gearboxFunction
    )
end
