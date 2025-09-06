mutable struct Gearbox
    ratio::carParameter
    torqueIn::carParameter
    angularFrequencyIn::carParameter
    torqueOut::carParameter
    angularFrequencyOut::carParameter
    loss::carParameter
    gearboxFunction
    
end


function createCTU25gearbox()
    ratio = carParameter(11.46,"Gearbox ratio","-")
    torqueIn = carParameter(0,"torque in","-")
    speedIn = carParameter(0,"velocity in","rad/s")

    torqueOut = carParameter(0,"torque out","-")
    speedOut = carParameter(0,"velocity out","rad/s")
    loss = carParameter(0,"loss","W")

    function gearboxFunction(gearbox) 
        gearbox.torqueOut = gearbox.torqueIn * 11.46
        gearbox.speedOut = gearbox.speedOut / 11/46
    end
    gearbox = Gearbox(
        ratio,
        torqueIn,
        speedIn,
        torqueOut,
        speedOut,
        loss,
        gearboxFunction
    )
end
