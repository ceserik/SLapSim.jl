function createLobattoIIIA_Adaptive(f, stages, model, nControls, nStates, track)
    tableau = TableauLobattoIIIA(stages)

    function createDynamicConstraints(segment_edges, initialization)
        number_of_segments = length(segment_edges) - 1
        totalPoints = number_of_segments * (stages - 1) + 1
        X = Matrix{VariableRef}(undef, totalPoints, nStates)
        U = Matrix{VariableRef}(undef, totalPoints, nControls)

        for i = 1:totalPoints
            for j = 1:nStates
                X[i, j] = @variable(model)
            end
            for j = 1:nControls
                U[i, j] = @variable(model)
            end
        end

        # Build full list s_all: each node followed by its collocation points (except last node)
        s_all = zeros(totalPoints)
        h_all = diff(segment_edges)
        idx = 1
        for segment = 1:number_of_segments #boundary nodes
            s_all[idx] = segment_edges[segment]
            h = segment_edges[segment+1] - segment_edges[segment]
            idx += 1

            for point = 2:(stages-1) #inner nodes
                #@infiltrate
                s_all[idx] = segment_edges[segment] + h * tableau.c[point]
                idx += 1
            end
        end

        s_all[end] = segment_edges[end]
        #create constraints using butcher tableau of lobattoIIIA
        segment_start_idx = 1
        f(X[segment_start_idx, :], U[segment_start_idx, :], s_all[segment_start_idx], model)
        for segment = 1:number_of_segments
            F = [f(X[segment_start_idx+col-1, :], U[segment_start_idx+col-1, :], s_all[segment_start_idx+col-1], nothing) for col in 1:stages]
            for stage = 2:stages
                fxSum = zeros(NonlinearExpr, nStates)
                for col = 1:stages
                    #@infiltrate     
                    fxSum += tableau.a[stage, col] * F[col]
                end
                ## quadrature constraint
                @constraint(model, X[segment_start_idx+stage-1, :] .== X[segment_start_idx, :] + h_all[segment] * fxSum)
                f(X[segment_start_idx+stage-1, :], U[segment_start_idx+stage-1, :], s_all[segment_start_idx+stage-1], model)
                #@infiltrate
            end
            segment_start_idx = segment_start_idx + stages - 1
        end

        #initialize all states and controls
        for idx = eachindex(s_all)
            init_vals = initialization.states(s_all[idx])
            init_u = initialization.controls(s_all[idx])
            set_start_value.(X[idx, :], init_vals)
            set_start_value.(U[idx, :], init_u)
        end

        return [model, X, U, s_all, segment_edges]
    end

    function create_interpolation(X_vals, U_vals, s_all, segment_edges)
        #shoul be real inteprolation but naaaaaaah, later
        function state_interp(s)
            node = searchsortedlast(s_all, s)
            return X_vals[node,:]
        end
        function control_interp(s)
            node = searchsortedlast(s_all, s)
            return U_vals[node,:]
        end
        return Result_interpolation(state_interp, control_interp, s_all)
    end

    LobattoIIIAMethod = Collocation(
        createDynamicConstraints,
        create_interpolation,
        nothing,
        nothing,
        f
    )
    return LobattoIIIAMethod
end