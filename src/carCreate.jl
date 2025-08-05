include("carExample.jl")
include("carParams.jl")

#Simple mass point model
function createCTU25()
    mass = carParameter(280.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    speed = carParameter(15.0,"Speed X","m/s")

    p = carParameters(
        mass,
        motorForce,
        CL,
        CD,
        speed,
    )

    function controlMapping(input, controls)
        input.motorForce.value = controls[1]  
        #paramcopy.CL.value = controls[2]
        return input
    end

    function stateMapping(input,states)
        input.speed.value = states[1]
        return input
    end
    car = Car(
        massPointCar,
        p,
        controlMapping,
        stateMapping,
    )

    return car
end


    
