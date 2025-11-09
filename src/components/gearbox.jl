
mutable struct Gearbox
    ratio::carParameter{carVar}
    torqueIn::carParameter{carVar}
    angularFrequencyIn::carParameter{carVar}
    torqueOut::carParameter{carVar}
    angularFrequencyOut::carParameter{carVar}
    loss::carParameter{carVar}
    f::Function
end

Base.show(io::IO, ::MIME"text/plain", obj::Gearbox) = prettyPrintComponent(io, obj)

#function Base.show(io::IO, ::MIME"text/plain", obj::Gearbox)
#    T = typeof(obj)
#    println(io, "$(T)(")
#    # Calculate maximum field name length for alignment
#    max_width = maximum(length(string(field)) for field in fieldnames(T))
#    for field in fieldnames(T)
#        # Pad field names with spaces for alignment
#        field_str = rpad(string(field), max_width)
#        println(io, "  $(field_str) = ", getfield(obj, field))
#    end
#    print(io, ")")
#end
#
function createCTU25gearbox()
    ratio = carParameter{carVar}(11.46,"Gearbox ratio","-")
    torqueIn = carParameter{carVar}(0.0,"torque in","-")
    speedIn = carParameter{carVar}(0.0,"velocity in","rad/s")

    torqueOut = carParameter{carVar}(0.0,"torque out","-")
    angularFrequencyOut = carParameter{carVar}(0.0,"velocity out","rad/s")
    loss = carParameter{carVar}(0.0,"loss","Watt")

    function gearboxFunction() 
        gearbox.torqueOut.value = gearbox.torqueIn.value * 11.46
        gearbox.angularFrequencyOut.value = gearbox.angularFrequencyOut.value / 11.46
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
