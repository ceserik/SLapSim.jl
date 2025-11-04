using SpecialPolynomials
using RungeKutta
using JuMP, Ipopt, Zygote
using Infiltrator
using HermiteInterpolation

struct Collocation
    createConstraints
    createInterpolator
    tableau
    interpolator::Union{Function,Nothing}
    f  
end


function createLobattoIIIA(stage,f)
    tableau = TableauLobattoIIIA(stage)
    stages = tableau.s
    function createDynamicConstraints(f, Xsize, Usize, iterpolationFunction, sampplingPoints, model, X_init, U_init)
        #create stages
        # potrebujem sampplingPoints + sampplingPoints* stage-2 dalsich
        #@infiltrate
        if (length(X_init) != length(U_init))
            U_init = [U_init; U_init[end, :]]
        end

        c = tableau.s
        stages = tableau.s
        N = length(sampplingPoints)

        numberOfCollocationPoints = N * (stages - 2)
        totalPoints = N + numberOfCollocationPoints

        @variable(model, X[i=1:totalPoints, j=1:Xsize]) #vyriesit inicializaciu
        @variable(model, U[i=1:totalPoints, j=1:Usize])

        #create vector of points which are node points (not collocation points)
        step4node = stages - 1
        Xnode = X[1:step4node:totalPoints, :]
        Unode = U[1:step4node:totalPoints, :]

        
        
        # Create correct sampling points
        s_all = zeros(totalPoints)
        k = 1
        for i = 1:N-1
            #@infiltrate
            h = iterpolationFunction(i+1) - iterpolationFunction(i)
            for stage = 1:stages -1
               
                s_all[k] = iterpolationFunction(i) + h*tableau.c[stage]
                #println(i)
                #println(stage)
                #println("")
                k +=1
            end
        end
        s_all[end] = iterpolationFunction(N)
        @infiltrate
        
        k = 1
        for xd = 1:N-1
            # Calculate the starting index for this interval
            interval_start_idx = 1 + (xd - 1) * (stages - 1)

            for stage = 2:tableau.s
                fxSum = zeros(NonlinearExpr, Xsize)
                for col = 1:tableau.s
                    # Correctly index into X for this interval
                    x_idx = interval_start_idx + col - 1
                    fxSum += tableau.a[stage, col] * f(X[x_idx, :], U[k, :],s_all[k])# sem pridat drahu
                end

                h = (iterpolationFunction(xd + 1) - iterpolationFunction(xd))

                # The next point in the sequence
                next_idx = interval_start_idx + stage - 1
                @constraint(model, X[next_idx, :] == X[interval_start_idx, :] + h * (fxSum))

                # set the intial solution for all points, even collocation points
                for i = 1:Xsize
                    s_stage = h * tableau.c[stage]
                    # @infiltrate
                    X_kInit = interp1(iterpolationFunction(sampplingPoints), X_init[:, i], iterpolationFunction(xd) + s_stage, "PCHIP")
                    set_start_value(X[next_idx, i], X_kInit)

                end

                for i = 1:Usize
                    s_stage = h * tableau.c[stage]
                    u_kInit = interp1(iterpolationFunction(sampplingPoints), U_init[:, i], iterpolationFunction(xd) + s_stage, "PCHIP")
                    set_start_value(U[next_idx, i], u_kInit)
                end
                k += 1
            end
        end
        return [model, X, U, Xnode, Unode[1:end-1, :],s_all] #this controls vector being shorter has to be better because track can be closed
    end


    function createLobattoInterpolator(x_values,u_values, s_all)
        step = stages - 1
        nodes_x = x_values[1:step:end, :]
        nodes_s = s_all[1:step:end, :]
        


        function interpolate(queryPoints)
            out = zeros(length(queryPoints),size(nodes_x,2))
            #identify segment  AKA node index
            for pointIDX = 1:length(queryPoints)
                segment_idx = -1
                #find correct segment
                for i = 1:length(nodes_s)
                    if nodes_s[i] > queryPoints[pointIDX]
                        segment_idx = i - 1
                        break
                    end
                end
                #@infiltrate
                segment_start = segment_idx * step -1
                x_segment = x_values[segment_start:segment_start+step, :]
                u_segment = u_values[segment_start:segment_start+step, :]
                s_segment = s_all[segment_start:segment_start+step]
                
                println(s_segment)
                println(x_segment)
                
                #if(length(s_segment) != length(x_segment))
                #    @infiltrate
                #end
                derivatives = zeros(stages,size(x_segment,2))
                for i = 1:stages
                    derivatives[i,:] =f(x_segment[i,:],u_segment[i,:],s_segment[i])
                end
                
                for x_idx = eachindex(x_segment[1,:])
                    polynom = HermiteInterpolation.fit(s_segment,x_segment[:,x_idx],derivatives[:,x_idx])
                    #@infiltrate
                    x_interp = polynom(queryPoints[pointIDX])
                    out[pointIDX,x_idx] = x_interp
                end
                #@infiltrate
                #poly = Lagrange(s_segment, vec(segment_values))
                #out[pointIDX] = poly(queryPoints[pointIDX])
                
            end
            #@infiltrate
            return out
        end

        return interpolate


    end


    LobattoIIIAMethod = Collocation(
        createDynamicConstraints,
        createLobattoInterpolator,
        tableau,
        nothing,
        f
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