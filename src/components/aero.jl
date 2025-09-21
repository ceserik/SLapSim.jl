mutable struct Aero
    CL::carParameter{carVar}
    CD::carParameter{carVar}
    CoP::carParameter{carVar}
end
#Base.show(io::IO, ::MIME"text/plain", obj) = prettyPrintComponent(io, obj)

function createBasicAero()
    CL = carParameter{carVar}(-5.0,"Lift coeffcient","-")
    CD = carParameter{carVar}(2.0,"Drag coeffcient","-")
    CoP = carParameter{carVar}(0.5,"Cener of pressure on front","-")
    aero = Aero(
        CL,
        CD,
        CoP
    )
    return aero

end
