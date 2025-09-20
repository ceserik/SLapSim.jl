mutable struct Aero
    CL::carParameter{Float64}
    CD::carParameter{Float64}
    CoP::carParameter{Float64}
end
Base.show(io::IO, ::MIME"text/plain", obj) = prettyPrintComponent(io, obj)

function createBasicAero()
    CL = carParameter(-5.0,"Lift coeffcient","-")
    CD = carParameter(2.0,"Drag coeffcient","-")
    CoP = carParameter(0.5,"Cener of pressure on front","-")
    aero = Aero(
        CL,
        CD,
        CoP
    )
    return aero

end
