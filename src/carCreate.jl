include("carExample.jl")
include("carParams.jl")

#Simple mass point model
function createCTU25()
    mass = carParameter(1500.0, "Mass", "kg")
    motorForce = carParameter(1000.0, "motorForce", "N")
    CL = carParameter(2.8, "Lift Coefficient", "-")
    CD = carParameter(1.0, "Drag Coefficient", "-")

    p = carParameters(
        mass,
        motorForce,
        CL,
        CD
    )

    function controlMapping(params, controls)
        paramcopy = deepcopy(params)
        paramcopy.motorForce.value = controls[1]  
        paramcopy.CL.value = controls[2]
        return paramcopy
    end

    car = Car(
        massPointCar,
        p,
        controlMapping
    )

    return car
end


    
