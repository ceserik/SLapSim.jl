using SpecialPolynomials
using RungeKutta
using JuMP, Ipopt, Zygote
using Infiltrator


struct Collocation
    createConstraints
    interpolate
    tableau
end


function createLobattoIIIA(stage)
    tableau = TableauLobattoIIIA(stage);
    function createDynamicConstraints(f,Xsize,Usize,independentVariable,sampplingPoints,model) 
        #create stages
        # potrebujem sampplingPoints + sampplingPoints* stage-2 dalsich
       

        c = tableau.s
        stages = tableau.s
        N = length(sampplingPoints)

        numberOfCollocationPoints = N*(stages - 2)

        @variable(model, X[i = 1:N + numberOfCollocationPoints,j = 1:Xsize]) #vyriesit inicializaciu
        @variable(model, U[i = 1:N + numberOfCollocationPoints,j = 1:Usize])


       
        for k = 1:N
            for stage = 2:tableau.s
                fxSum = carVar
                #hStage = tableau.c[stage]*h
                
                for col = 1:tableau.s
                    @infiltrate
                    fxSum += tableau.a[stage,col] * f(X[k,:],U[k,:]) #construct etimation of collocation point
                end
                
                h = sampplingPoints[k+1] - sampplingPoints[k]
                nextPoint = k*(stages - 2) +stage +1
                @constraint(model, X[k+1+stage,:]== X[k,:] + f(X[k,:] + h*(fxSum),U[k,:]) ) 

            end
        end
        return model
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

stage = 2
   function f(z, u)
        return [u; z[1]; cos(z[2]); sin(z[2])]
    end
model = JuMP.Model(Ipopt.Optimizer)

xdlol = createLobattoIIIA(stage)

xdlol.createConstraints(f,2,1,0,track.sampleDistances,model)