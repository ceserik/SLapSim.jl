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
    #tableau = TableauRunge()
    stages = tableau.s
    function createDynamicConstraints(f, Xsize::Int64, Usize::Int64, iterpolationFunction, sampplingPoints::Vector{Float64}, model::JuMP.Model, X_init, U_init)
        # Ensure U_init has same number of node entries as X_init (extend last row if necessary)
        if size(U_init, 1) != size(X_init, 1)

            U_init = vcat(U_init, U_init[end, :]')
        end

        stages = tableau.s
        N = length(sampplingPoints)

        numberOfCollocationPoints = (N - 1) * (stages - 2)
        totalPoints = N + numberOfCollocationPoints

        # Create variables sequentially (one time point at a time) for proper Jacobian ordering
        # This ensures lower triangular structure in the Jacobian
        X = Matrix{VariableRef}(undef, totalPoints, Xsize)
        U = Matrix{VariableRef}(undef, totalPoints, Usize)
        
        for i = 1:totalPoints
            for j = 1:Xsize
                X[i, j] = @variable(model)
            end
            for j = 1:Usize
                U[i, j] = @variable(model)
            end
        end

        # node indices stride
        step4node = stages - 1
        Xnode = X[1:step4node:totalPoints, :]
        Unode = U[1:step4node:totalPoints, :]

        # Build node sampling positions (apply iterpolationFunction to sampplingPoints once)
        node_s = sampplingPoints
        # Build full list s_all: each node followed by its collocation points (except last node)
        s_all = zeros(totalPoints)
        idx = 1
        for i = 1:N
            s_all[idx] = node_s[i]
            idx += 1
            if i < N
                h_node = node_s[i + 1] - node_s[i]
                for st = 2:(stages - 1)
                    s_all[idx] = node_s[i] + h_node * tableau.c[st]
                    idx += 1
                end
            end
        end
        constrained = falses(totalPoints)
        # Collocation / continuity constraints for each interval
        for interval = 1:(N - 1)
            interval_start_idx = 1 + (interval - 1) * (stages - 1)
            node_start_s = node_s[interval]
            node_end_s = node_s[interval + 1]
            h = node_end_s - node_start_s

            for stage = 2:stages
                # index in X/U corresponding to this collocation/stage point
                next_idx = interval_start_idx + stage - 1

                
                # Build Runge-Kutta weighted sum of f evaluated at the stage points in this interval
                fxSum = zeros(NonlinearExpr, Xsize)
                for col = 1:stages
                    x_idx = interval_start_idx + col - 1
                    fxSum += tableau.a[stage, col] * f(X[x_idx, :], U[x_idx, :], s_all[x_idx],false)
                    #xdddd = f(X[x_idx, :], U[x_idx, :], s_all[x_idx],true)
                end

                 # Only add model constraints for each time point once
                if !constrained[next_idx]
                    _ = f(X[next_idx, :], U[next_idx, :], s_all[next_idx], true)
                    constrained[next_idx] = true
                end
                if !constrained[interval_start_idx]
                    _ = f(X[interval_start_idx, :], U[interval_start_idx, :], s_all[interval_start_idx], true)
                    constrained[interval_start_idx] = true
                end

                ## quadrature constraint
                @constraint(model, X[next_idx, :] .== X[interval_start_idx, :] + h * fxSum)

                # set the initial guess (start value) for this collocation point using interpolation over node initial guesses
                s_stage = node_start_s + h * tableau.c[stage]
                for i = 1:Xsize
                    X_kInit = interp1(node_s, X_init[:, i], s_stage, "PCHIP")
                    set_start_value(X[next_idx, i], X_kInit)
                    #println(X_kInit)

                end
                #println()
                for i = 1:Usize
                    u_kInit = interp1(node_s, U_init[:, i], s_stage, "PCHIP")
                    set_start_value(U[next_idx, i], u_kInit)
                    #println(u_kInit)
                end
            end
        end
        
        for i = 1:Xsize
            set_start_value(X[1, i], X_init[1, i])
        end

        for i = 1:Usize
            set_start_value(U[1, i], U_init[1, i])
        end

        return [model, X, U, Xnode, Unode[1:end-1, :], s_all]
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

                # Handle edge cases
                if segment_idx == -1
                    if queryPoints[pointIDX] < nodes_s[1]
                        segment_idx = 1
                    else
                        segment_idx = length(nodes_s) - 1
                    end
                end
                if segment_idx ==0
                    segment_idx =1
                end
                segment_start = 1 + (segment_idx - 1) * step
                x_segment = x_values[segment_start:segment_start+step, :]
                u_segment = u_values[segment_start:segment_start+step, :]
                s_segment = s_all[segment_start:segment_start+step]
                
                derivatives = zeros(stages,size(x_segment,2))
                for i = 1:stages
                    derivatives[i,:] =f(x_segment[i,:],u_segment[i,:],s_segment[i],false)
                end
                
                for x_idx = eachindex(x_segment[1,:])
                    polynom = HermiteInterpolation.fit(s_segment,x_segment[:,x_idx],derivatives[:,x_idx])
                    x_interp = polynom(queryPoints[pointIDX])
                    out[pointIDX,x_idx] = x_interp
                end
                
            end
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
