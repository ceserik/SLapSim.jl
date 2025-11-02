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
    function createDynamicConstraints(f,Xsize,Usize,independentVariable,sampplingPoints,model,X_init,Y_init) 
        #create stages
        # potrebujem sampplingPoints + sampplingPoints* stage-2 dalsich
       
        c = tableau.s
        stages = tableau.s
        N = length(sampplingPoints)

        numberOfCollocationPoints = N*(stages - 2)
        totalPoints = N + numberOfCollocationPoints
        

        
        @variable(model, X[i = 1:totalPoints,j = 1:Xsize]) #vyriesit inicializaciu
        @variable(model, U[i = 1:totalPoints,j = 1:Usize])

        #create vector of points which are node points (not collocation points)
        step4node = stages-1
        Xnode = X[1:step4node:totalPoints,:]
        Unode = U[1:step4node:totalPoints,:]

        k = 1
        for xd = 1:N-1
            # Calculate the starting index for this interval
            interval_start_idx = 1 + (xd-1)*(stages-1)
            
            for stage = 2:tableau.s
                fxSum = zeros(NonlinearExpr, Xsize)
                for col = 1:tableau.s
                    # Correctly index into X for this interval
                    x_idx = interval_start_idx + col - 1
                    fxSum += tableau.a[stage,col] * f(X[x_idx,:],U[k,:])
                end
                
                h = (sampplingPoints[xd+1] - sampplingPoints[xd])
                
                # The next point in the sequence
                next_idx = interval_start_idx + stage - 1
                @constraint(model, X[next_idx,:] == X[interval_start_idx,:] + h*(fxSum))

                for i = 1:Xsize
                    h_stage = h * tableau.c[stage]
                    X_kInit = interp1(sampplingPoints, X_init[:,i], sampplingPoints[xd] + h_stage, "PCHIP")
                    set_start_value(X[next_idx,i], X_kInit)
                end
                k += 1
            end
        end
        return  [model,  X, U[1:end-1,:], Xnode,Unode[1:end-1,:]] #this controls vector being shorter has to be better because track can be closed
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

#stage = 5
#   function f(z, u)
#        return [u; z[1]; cos(z[2]); sin(z[2])]
#    end
#model = JuMP.Model(Ipopt.Optimizer)
#
#xdlol = createLobattoIIIA(stage)
#
#xdlol.createConstraints(f,4,1,0,track.sampleDistances,model)