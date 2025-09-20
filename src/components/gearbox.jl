
mutable struct Gearbox
    ratio::carParameter{Float64}
    torqueIn::carParameter{Float64}
    angularFrequencyIn::carParameter{Float64}
    torqueOut::carParameter{Float64}
    angularFrequencyOut::carParameter{Float64}
    loss::carParameter{Float64}
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
    ratio = carParameter(11.46,"Gearbox ratio","-")
    torqueIn = carParameter(0.0,"torque in","-")
    speedIn = carParameter(0.0,"velocity in","rad/s")

    torqueOut = carParameter(0.0,"torque out","-")
    angularFrequencyOut = carParameter(0.0,"velocity out","rad/s")
    loss = carParameter(0.0,"loss","Watt")

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
