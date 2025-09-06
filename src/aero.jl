mutable struct Aero
    CL::carParameter
    CD::carParameter
    CoP::carParameter
end


function createBasicAero()
    CL = carParameter(-5.0,"Lift coeffcient","-")
    CD = carParameter(2,"Drag coeffcient","-")
    CoP = carParameter(0.5,"Cener of pressure on front","-")
    aero = Aero(
        CL,
        CD,
        CoP
    )
    return aero

end
