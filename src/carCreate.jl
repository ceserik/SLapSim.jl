include("carExample.jl")
include("carParams.jl")

#Simple mass point model
function createCTU25()
    mass = carParameter(150.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    CL = carParameter(2.8, "Lift Coefficient", "-")
    CD = carParameter(1.0, "Drag Coefficient", "-")
    speed = carParameter(15.0,"Speed X","m/s")

    p = carParameters(
        mass,
        motorForce,
        CL,
        CD,
        speed,
    )

    function controlMapping(params, controls)
        paramcopy = deepcopy(params)
        paramcopy.motorForce.value = controls[1]  
        #paramcopy.CL.value = controls[2]
        return paramcopy
    end

    function stateMapping(params,states)
        paramcopy = deepcopy(params)
        paramcopy.speed.value = states[1]
        return paramcopy
    end
    car = Car(
        massPointCar,
        p,
        controlMapping,
        stateMapping,
    )

    return car
end


    
