function createLobattoIIIA_Adaptive(f, stages, model, nControls, nStates, track;
                                    x_scale=ones(nStates), u_scale=ones(nControls),
                                    x_lb=-x_scale, x_ub=x_scale,
                                    u_lb=-u_scale, u_ub=u_scale)
    tableau = TableauLobattoIIIA(stages)
    τ_ref  = tableau.c
    w_bary = barycentric_weights(τ_ref)

    function createDynamicConstraints(segment_edges, initialization)
        number_of_segments = length(segment_edges) - 1
        totalPoints = number_of_segments * (stages - 1) + 1

        X_raw = Matrix{VariableRef}(undef, totalPoints, nStates)
        U_raw = Matrix{VariableRef}(undef, totalPoints, nControls)
        x_lb_raw = x_lb ./ x_scale
        x_ub_raw = x_ub ./ x_scale
        u_lb_raw = u_lb ./ u_scale
        u_ub_raw = u_ub ./ u_scale
        println("State bounds & scales:")
        for j in eachindex(x_lb)
            println("  x[$j]: lb=$(round(x_lb[j],digits=3))  ub=$(round(x_ub[j],digits=3))  scale=$(round(x_scale[j],digits=3))")
        end
        println("Control bounds & scales:")
        for j in eachindex(u_lb)
            println("  u[$j]: lb=$(round(u_lb[j],digits=3))  ub=$(round(u_ub[j],digits=3))  scale=$(round(u_scale[j],digits=3))")
        end

        for i = 1:totalPoints
            for j = 1:nStates
                X_raw[i, j] = @variable(model, lower_bound=x_lb_raw[j], upper_bound=x_ub_raw[j])
            end
            for j = 1:nControls
                U_raw[i, j] = @variable(model, lower_bound=u_lb_raw[j], upper_bound=u_ub_raw[j])
            end
        end

        #scaling
        X = X_raw .* x_scale'
        U = U_raw .* u_scale'

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

        #initialize all states and controls (set on raw O(1) variables)
        for idx = eachindex(s_all)
            init_vals = initialization.states(s_all[idx])
            init_u = initialization.controls(s_all[idx])
            set_start_value.(X_raw[idx, :], init_vals ./ x_scale)
            set_start_value.(U_raw[idx, :], init_u ./ u_scale)
        end

        return [model, X, U, s_all, segment_edges]
    end
 
    function create_interpolation(X_vals, U_vals, s_all, segment_edges)
        number_of_segments = length(segment_edges) - 1
        stride = stages - 1  # points per segment (adjacent segments share boundary)
 
        function find_segment(s)
            for i in 1:number_of_segments
                s <= segment_edges[i+1] && return i
            end
            return number_of_segments
        end
 

        function state_interp(s)
            seg   = find_segment(s)
            h     = segment_edges[seg+1] - segment_edges[seg]
            τ_eval = (s - segment_edges[seg]) / h       # map to [0,1]
            i0    = (seg - 1) * stride + 1
            i_end = i0 + stages - 1
            x_seg = X_vals[i0:i_end, :]
            return [_bary_eval(τ_ref, w_bary, x_seg[:, k], τ_eval) for k in 1:nStates]
        end
 
        function control_interp(s)
            seg    = find_segment(s)
            h      = segment_edges[seg+1] - segment_edges[seg]
            τ_eval = (s - segment_edges[seg]) / h
            i0     = (seg - 1) * stride + 1
            i_end  = i0 + stages - 1
            u_seg  = U_vals[i0:i_end, :]
            return [_bary_eval(τ_ref, w_bary, u_seg[:, k], τ_eval) for k in 1:nControls]
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