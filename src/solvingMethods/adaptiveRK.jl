function createLobattoIIIA_Adaptive(f, stages, model, nControls, nStates, track)
    tableau = TableauLobattoIIIA(stages)
    τ_ref  = tableau.c
    w_bary = barycentric_weights(τ_ref)

    function createDynamicConstraints(segment_edges, initialization)
        number_of_segments = length(segment_edges) - 1
        totalPoints = number_of_segments * (stages - 1) + 1

        # Build s_all FIRST so per-node bounds can be evaluated during variable creation.
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

        X = Matrix{VariableRef}(undef, totalPoints, nStates)
        U = Matrix{VariableRef}(undef, totalPoints, nControls)

        for i = 1:totalPoints
            for j = 1:nStates
                X[i, j] = @variable(model, lower_bound=-1.0, upper_bound=1.0)
            end
            for j = 1:nControls
                U[i, j] = @variable(model, lower_bound=-1.0, upper_bound=1.0)
            end
        end

        # === Intermediate-F variant (comment this block out to revert to inline f) ===
        F_raw = Matrix{VariableRef}(undef, totalPoints, nStates)
        for i = 1:totalPoints, j = 1:nStates
            F_raw[i, j] = @variable(model)
        end
        for i = 1:totalPoints
            dx_i = f(X[i, :], U[i, :], s_all[i], model)
            @constraint(model, F_raw[i, :] .== dx_i)
        end
        segment_start_idx = 1
        for segment = 1:number_of_segments
            for stage = 2:stages
                fxSum = @variable(model, [1:nStates])
                @constraint(model, fxSum .== sum(tableau.a[stage, col] .* F_raw[segment_start_idx+col-1, :] for col in 1:stages))

                @constraint(model, X[segment_start_idx+stage-1, :] .== X[segment_start_idx, :] .+ h_all[segment] .* fxSum)
            end
            segment_start_idx = segment_start_idx + stages - 1
        end
        # === end intermediate-F variant ===

        # Original inline variant (kept for easy toggle):
        #segment_start_idx = 1
        #f(X[segment_start_idx, :], U[segment_start_idx, :], s_all[segment_start_idx], model)
        #for segment = 1:number_of_segments
        #    F = [f(X[segment_start_idx+col-1, :], U[segment_start_idx+col-1, :], s_all[segment_start_idx+col-1], nothing) for col in 1:stages]
        #    # === odow intermediate-F toggle: uncomment 2 lines below to split dynamics into aux vars (sparser hessian, ~2× faster) ===
        #    #F_aux = [@variable(model) for _=1:nStates, _=1:stages]
        #    #@constraint(model, F_aux .== hcat(F...)); F = [F_aux[:,col] for col in 1:stages]
        #    for stage = 2:stages
        #        fxSum = zeros(NonlinearExpr, nStates)
        #        for col = 1:stages
        #            fxSum += tableau.a[stage, col] * F[col]
        #        end
        #        @constraint(model, X_raw[segment_start_idx+stage-1, :] .== X_raw[segment_start_idx, :] .+ h_all[segment] .* fxSum ./ x_scale)
        #        f(X[segment_start_idx+stage-1, :], U[segment_start_idx+stage-1, :], s_all[segment_start_idx+stage-1], model)
        #    end
        #    segment_start_idx = segment_start_idx + stages - 1
        #end

        #initialize all states and controls 
        for idx = eachindex(s_all)
            init_vals = initialization.states(s_all[idx])
            init_u = initialization.controls(s_all[idx])
            set_start_value.(X[idx, :], init_vals)
            set_start_value.(U[idx, :], init_u)
        end

        return [model, X, U, s_all, segment_edges, s_all]
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
 
        return Result_interpolation(state_interp, control_interp, s_all, collect(Float64, segment_edges))
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