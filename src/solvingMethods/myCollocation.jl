using SpecialPolynomials
using RungeKutta
using JuMP
using Infiltrator
struct Collocation
    createConstraints
    interpolate
    tableau

end




function createLobattoIIIA(stage)
    tableau = TableauLobattoIIIA(stage)
    function createDynamicConstraints(f,X,U,model)
        #create stages
        N = len(X)
        A = zeros(tableau.s,1)
        h = s[k+1] - s[k]
        c = tableau.s
        #@variable(model,Xc[1:N,1:tab.s]) # Xc ako collocation

        for stage = 2:tableau.s
            fxSum = carVar
            hStage = tableau.c[stage]*h
            for col = 1:tableau.s
                fxSum += tableau.a[stage,col] * f(X[k,:],U[k,:])
            end
            
            @constraint(model, X[k+1,:]== X[k+1,:] + f(X[k,:] + h*(fxSum),U[k,:]) ) 

        end

    end

    function LobattoInterpolation(points,querypoints)
        #@infiltrate
        p =  Lagrange(collect(tableau.c),points)
        p.(querypoints)
    end


    LobattoIIIAMethod = Collocation(
        createDynamicConstraints,
        LobattoInterpolation,
        tableau
    )
    return LobattoIIIAMethod
end