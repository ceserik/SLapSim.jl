using SpecialPolynomials
using RungeKutta
using JuMP, Ipopt, Zygote
using Infiltrator


struct Collocation
    createConstraints
    createInterpolator
    tableau
    interpolator::Union{Function, Nothing}  # Store the actual interpolator here
end


function createLobattoIIIA(stage)
    tableau = TableauLobattoIIIA(stage);
    stages = tableau.s  
    function createDynamicConstraints(f,Xsize,Usize,iterpolationFunction,sampplingPoints,model,X_init,U_init) 
        #create stages
        # potrebujem sampplingPoints + sampplingPoints* stage-2 dalsich
        #@infiltrate
        if (length(X_init) != length(U_init))
            U_init = [U_init; U_init[end,:] ]
        end


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
                    fxSum += tableau.a[stage,col] * f(X[x_idx,:],U[k,:])# sem pridat drahu
                end
                
                h = (iterpolationFunction(xd+1) - iterpolationFunction(xd))
                
                # The next point in the sequence
                next_idx = interval_start_idx + stage - 1
                @constraint(model, X[next_idx,:] == X[interval_start_idx,:] + h*(fxSum))
                
                # set the intial solution for all points, even collocation points
                for i = 1:Xsize
                    s_stage = h * tableau.c[stage]
                   # @infiltrate
                    X_kInit = interp1(iterpolationFunction(sampplingPoints), X_init[:,i], iterpolationFunction(xd) + s_stage, "PCHIP")
                    set_start_value(X[next_idx,i], X_kInit)
                    
                end

                for i = 1:Usize
                    s_stage = h * tableau.c[stage]
                    u_kInit = interp1(iterpolationFunction(sampplingPoints), U_init[:,i], iterpolationFunction(xd) + s_stage, "PCHIP")
                    set_start_value(U[next_idx,i], u_kInit)
                end
                k += 1
            end
        end
        return  [model,  X, U, Xnode,Unode[1:end-1,:]] #this controls vector being shorter has to be better because track can be closed
    end


        function createLobattoInterpolator(X_values, samplingDistances)
        # Return a function that only needs querypoints
        function interpolate(querypoints)
            nIntervals = length(samplingDistances) - 1
            nStates = size(X_values, 2)
            result = zeros(length(querypoints), nStates)
            
            for (idx, s_query) in enumerate(querypoints)
                # Find interval
                interval = clamp(searchsortedlast(samplingDistances, s_query), 1, nIntervals)
                
                # Local coordinate in [0, 1]
                s_start = samplingDistances[interval]
                s_end = samplingDistances[interval + 1]
                tau = (s_query - s_start) / (s_end - s_start)
                
                # Extract collocation points for this interval
                interval_start = 1 + (interval - 1) * (stages - 1)
                points = X_values[interval_start:(interval_start + stages - 1), :]
                
                # Lagrange interpolation at tau
                for state = 1:nStates
                    poly = Lagrange(collect(tableau.c), points[:, state])
                    result[idx, state] = poly(tau)
                end
            end
            
            return result
        end
        
        return interpolate
    end


    LobattoIIIAMethod = Collocation(
        createDynamicConstraints,
        createLobattoInterpolator,
        tableau,
        nothing
    )
    return LobattoIIIAMethod
end



#Hermite polynomial

#stage = 5
#   function f(z, u)
#        return [u; z[1]; cos(z[2]); sin(z[2])]
#    end
#model = JuMP.Model(Ipopt.Optimizer)
#
#xdlol = createLobattoIIIA(stage)
#
#xdlol.createConstraints(f,4,1,0,track.sampleDistances,model)