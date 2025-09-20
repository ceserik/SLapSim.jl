mutable struct Chassis
    mass::carParameter{Float64}
    hitbox::Function
end

function createCTU25chassis()
    mass = carParameter(180.0,"mass","kg")

    function hitbox(car::Car,track::Track,model::JuMP.Model)
        @constraint(model,car.n .>= -1.5)
        @constraint(model,car.n .<= 1.5)
    end

    chassis = Chassis(
        mass,
        hitbox
    )
    return chassis
end
