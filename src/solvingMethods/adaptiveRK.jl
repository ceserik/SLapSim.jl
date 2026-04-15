function createLobattoIIIA_Adaptive(f, stages, model, nControls, nStates, track;
                                    x_scale=ones(nStates), u_scale=ones(nControls),
                                    x_lb=s -> -x_scale, x_ub=s -> x_scale,
                                    u_lb=s -> -u_scale, u_ub=s -> u_scale)
    tableau = TableauLobattoIIIA(stages)
    τ_ref  = tableau.c
    w_bary = barycentric_weights(τ_ref)

    # Evaluate a bound function at path position and verify that it returns the expected number of entries.
    function _eval_bound(bound_fn, s, expected_length)
        bound_values = bound_fn(s)
        length(bound_values) == expected_length || error(
            "bound function returned vector of length $(length(bound_values)), " *
            "expected $(expected_length) at s=$(s)")
        return bound_values
    end

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

        X_raw = Matrix{VariableRef}(undef, totalPoints, nStates)
        U_raw = Matrix{VariableRef}(undef, totalPoints, nControls)

        # Print bounds & scales summary using the value at s=0 (representative).
        let s0 = s_all[1]
            x_lb0 = _eval_bound(x_lb, s0, nStates)
            x_ub0 = _eval_bound(x_ub, s0, nStates)
            u_lb0 = _eval_bound(u_lb, s0, nControls)
            u_ub0 = _eval_bound(u_ub, s0, nControls)
            println("State bounds & scales (at s=$(round(s0,digits=3))):")
            for j in eachindex(x_lb0)
                println("  x[$j]: lb=$(round(x_lb0[j],digits=3))  ub=$(round(x_ub0[j],digits=3))  scale=$(round(x_scale[j],digits=3))")
            end
            println("Control bounds & scales (at s=$(round(s0,digits=3))):")
            for j in eachindex(u_lb0)
                println("  u[$j]: lb=$(round(u_lb0[j],digits=3))  ub=$(round(u_ub0[j],digits=3))  scale=$(round(u_scale[j],digits=3))")
            end
        end

        for i = 1:totalPoints
            si = s_all[i]
            x_lb_i = _eval_bound(x_lb, si, nStates) ./ x_scale
            x_ub_i = _eval_bound(x_ub, si, nStates) ./ x_scale
            u_lb_i = _eval_bound(u_lb, si, nControls) ./ u_scale
            u_ub_i = _eval_bound(u_ub, si, nControls) ./ u_scale
            for j = 1:nStates
                X_raw[i, j] = @variable(model, lower_bound=x_lb_i[j], upper_bound=x_ub_i[j])
            end
            for j = 1:nControls
                U_raw[i, j] = @variable(model, lower_bound=u_lb_i[j], upper_bound=u_ub_i[j])
            end
        end

        #scaling
        X = X_raw .* x_scale'
        U = U_raw .* u_scale'

        # === Intermediate-F variant (comment this block out to revert to inline f) ===
        # Split dynamics: F_raw[i,:] == f(X,U,s)/x_scale as equality, then quadrature
        # is linear in F_raw → sparse hessian rows, faster AD (odow tip, ExaModels pattern).
        #F_raw = Matrix{VariableRef}(undef, totalPoints, nStates)
        #for i = 1:totalPoints, j = 1:nStates
        #    F_raw[i, j] = @variable(model)
        #end
        #for i = 1:totalPoints
        #    dx_i = f(X[i, :], U[i, :], s_all[i], model)
        #    @constraint(model, F_raw[i, :] .== dx_i ./ x_scale)
        #end
        #segment_start_idx = 1
        #for segment = 1:number_of_segments
        #    for stage = 2:stages
        #        fxSum = sum(tableau.a[stage, col] .* F_raw[segment_start_idx+col-1, :] for col in 1:stages)
        #        @constraint(model, X_raw[segment_start_idx+stage-1, :] .== X_raw[segment_start_idx, :] .+ h_all[segment] .* fxSum)
        #    end
        #    segment_start_idx = segment_start_idx + stages - 1
        #end
        # === end intermediate-F variant ===

        # Original inline variant (kept for easy toggle):
        segment_start_idx = 1
        f(X[segment_start_idx, :], U[segment_start_idx, :], s_all[segment_start_idx], model)
        for segment = 1:number_of_segments
            F = [f(X[segment_start_idx+col-1, :], U[segment_start_idx+col-1, :], s_all[segment_start_idx+col-1], nothing) for col in 1:stages]
            for stage = 2:stages
                fxSum = zeros(NonlinearExpr, nStates)
                for col = 1:stages
                    fxSum += tableau.a[stage, col] * F[col]
                end
                @constraint(model, X_raw[segment_start_idx+stage-1, :] .== X_raw[segment_start_idx, :] .+ h_all[segment] .* fxSum ./ x_scale)
                f(X[segment_start_idx+stage-1, :], U[segment_start_idx+stage-1, :], s_all[segment_start_idx+stage-1], model)
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