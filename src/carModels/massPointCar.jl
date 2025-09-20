
function massPointCar(car::Car,track::Track,k::Int64, optiModel::Union{Nothing, JuMP.Model} = nothing)
    # Simple  Mass point friction circle
    # inputs is a struct of all posible inputs of vehicle, controls are not separated from inputs
    # Get inputs
    m          = car.carParameters.mass.value
    inputForce = car.carParameters.motorForce.value
    CL         = car.carParameters.CL.value
    CD         = car.carParameters.CD.value
    maxPower   = car.carParameters.powerLimit.value

    # Get state
    vx         = car.carParameters.velocity.value[1]
        # this simplificqation is necceseary so that  universal function time -> path can be used
        #it means that car is always pointing the sma edirection as track and no sideways motion is possible
    car.carParameters.psi.value = track.theta[k]
    
    # Get track inputs
    c          = track.curvature[k]
    rho        = track.rho[k]
    μ          = track.μ[k]
    # Calculate forces

    Fz = 1/2 * rho * CL * vx^2
    Fy = m*vx^2*c

    maxMotorForce = 6507 #calculated for ctu25 should be added to inputs, vx torque characerisitic
    FxPowerMax = maxPower/(vx)
    #print(Fy)
    FxMaxsquared = max((Fz*μ + m*9.81*μ)^2 - Fy^2,0)
    # Add optimization constraints if model is provided

    if optiModel !== nothing
        ## tu urobit funkciu do ktorej dam obmedzenie a ona mi o spravi aby som nemusel stale davat ify
        ## aj ked to mozno je jedno lebo tento mass point bude mozno malo pouzivany?
        ## ale radsej to spravit, nech netreba prepisovat model lebo sa z toho zblaznim ked tam bude nieco inak
        @constraint(optiModel, (Fy/maxMotorForce)^2 + (inputForce/maxMotorForce)^2 <= ((Fz*μ + m*9.81*μ)/maxMotorForce)^2)
        @constraint(optiModel,inputForce/FxPowerMax <= FxPowerMax/FxPowerMax)
    else
        inputForce = min(inputForce, sqrt(FxMaxsquared))
        inputForce = max(inputForce, -sqrt(FxMaxsquared))
        inputForce = min(inputForce,maxMotorForce)
        inputForce = min(inputForce,FxPowerMax)
    end         

   # lessContraint((Fy/maxMotorForce)^2 + (inputForce/maxMotorForce)^2 ,((Fz*μ + m*9.81*μ)/maxMotorForce)^2,optiModel)
   # lessContraint(inputForce/FxPowerMax , FxPowerMax/FxPowerMax,optiModel)

        dstates =  [(inputForce -(1/2 *rho*CD*vx^2))/m 0  0    0     0 0]
    return dstates
end


#Simple mass point model
function createCTU25_1D()
    mass = carParameter(280.0, "Mass", "kg")
    inertia = carParameter(100.0, "Inertia", "kg*m^2")
    motorForce = carParameter(1000.0, "motorForce", "N")
    lateralForce = carParameter(0.0, "lateral Force", "N")
    CL = carParameter(5.0, "Lift Coefficient", "-")
    CD = carParameter(2.0, "Drag Coefficient", "-")
    velocity = carParameter([15.0,0,0],"Speed X","m/s")
    angularVelocity = carParameter([0.0, 0.0, 0.0], "angular velocity", "rad/s");
    powerLimit = carParameter(80000.0,"PowerLimit","W")
    vy = carParameter(0.0,"Speed Y","m/s")
    psi = carParameter(0.0,"heading","rad")
    n = carParameter(0.0,"Distance from centerline","m")
    nControls = carParameter(1.0,"number of controlled parameters","-")
    nStates = carParameter(1.0,"number of car states","-")

    p = carParameters(
        mass,
        inertia,
        motorForce,
        CL,
        CD,
        velocity,
        angularVelocity,
        psi,
        n,
        powerLimit,
        lateralForce,
        nControls,
        nStates
    )

    
    function controlMapping(car::Car,u::Union{Vector{VariableRef},Vector{Float64}})
        v = u[1]
        if isa(v, Number)
            # update existing numeric parameter value in-place
            car.carParameters.motorForce.value = float(v)
        else
            # replace parameter with one that stores the JuMP variable/expression
            car.carParameters.motorForce = carParameter(v, "motorForce", "N")
        end
        # other mappings if needed...
    end

    function mapping(car::Car,u::Union{Vector{VariableRef},Vector{Float64}},x::Union{Vector{VariableRef},Vector{Float64}})
        car.carParameters.motorForce.value = u[1]
        car.carParameters.vx.value = x[1]
    end

    function stateMapping(car::Car,x::Union{Vector{VariableRef},Vector{Float64}})
        # Build a 3-element vector of state components (pad with 0.0 if missing)
        vals = (
            length(x) >= 1 ? x[1] : 0.0,
            length(x) >= 2 ? x[2] : 0.0,
            length(x) >= 3 ? x[3] : 0.0
        )

        curp = car.carParameters.velocity

        # If all state entries are numeric and current storage is a numeric vector, update in-place
        if all(isa.(Tuple(vals), Number)) && isa(curp.value, AbstractVector)
            @inbounds begin
                n = min(length(curp.value), 3)
                for i in 1:n
                    curp.value[i] = float(vals[i])
                end
            end
        else
            # Otherwise replace the carParameter with one that holds the (possibly non-numeric) vector
            # Use a 1D vector constructor with commas so we don't create a 1×3 Matrix
            newvec = [vals[1], vals[2], vals[3]]
            car.carParameters.velocity = carParameter(newvec, curp.name, curp.unit)
        end
    end
    car = Car(
        massPointCar,
        p,
        controlMapping,
        stateMapping,
        mapping
    )

    return car
end

