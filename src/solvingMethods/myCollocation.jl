using Revise
using RungeKutta
using JuMP, Ipopt, Zygote
using Infiltrator
using HermiteInterpolation
using FastGaussQuadrature

struct Collocation
    createConstraints
    createInterpolator
    tableau
    interpolator::Union{Function,Nothing}
    f
end

# This function creates Legendre-Gauss-Lobatto Collocation
function createLobattoIIIA(stage, f)
    tableau = TableauLobattoIIIA(stage)
    stages = tableau.s
    function createDynamicConstraints(f, Xsize::Int64, Usize::Int64, iterpolationFunction, sampplingPoints::Vector{Float64}, model::JuMP.Model, X_init, U_init)
        if size(U_init, 1) != size(X_init, 1)
            U_init = vcat(U_init, U_init[end, :]')
        end

        stages = tableau.s
        N = length(sampplingPoints)

        numberOfCollocationPoints = (N - 1) * (stages - 2)
        totalPoints = N + numberOfCollocationPoints

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

        step4node = stages - 1
        Xnode = X[1:step4node:totalPoints, :]
        Unode = U[1:step4node:totalPoints, :]

        node_s = sampplingPoints
        s_all = zeros(totalPoints)
        idx = 1
        for i = 1:N
            s_all[idx] = node_s[i]
            idx += 1
            if i < N
                h_node = node_s[i+1] - node_s[i]
                for st = 2:(stages-1)
                    s_all[idx] = node_s[i] + h_node * tableau.c[st]
                    idx += 1
                end
            end
        end
        constrained = falses(totalPoints)
        for interval = 1:(N-1)
            interval_start_idx = 1 + (interval - 1) * (stages - 1)
            node_start_s = node_s[interval]
            node_end_s = node_s[interval+1]
            h = node_end_s - node_start_s

            for stage = 2:stages
                next_idx = interval_start_idx + stage - 1

                fxSum = zeros(NonlinearExpr, Xsize)
                for col = 1:stages
                    x_idx = interval_start_idx + col - 1
                    fxSum += tableau.a[stage, col] * f(X[x_idx, :], U[x_idx, :], s_all[x_idx], nothing)
                end

                if !constrained[next_idx]
                    _ = f(X[next_idx, :], U[next_idx, :], s_all[next_idx], model)
                    constrained[next_idx] = true
                end
                if !constrained[interval_start_idx]
                    _ = f(X[interval_start_idx, :], U[interval_start_idx, :], s_all[interval_start_idx], model)
                    constrained[interval_start_idx] = true
                end

                @constraint(model, X[next_idx, :] .== X[interval_start_idx, :] + h * fxSum)

                s_stage = node_start_s + h * tableau.c[stage]
                for i = 1:Xsize
                    X_kInit = interp1(node_s, X_init[:, i], s_stage, "PCHIP")
                    set_start_value(X[next_idx, i], X_kInit)
                end
                for i = 1:Usize
                    u_kInit = interp1(node_s, U_init[:, i], s_stage, "PCHIP")
                    set_start_value(U[next_idx, i], u_kInit)
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

    function createLobattoInterpolator(x_values, u_values, s_all)
        step = stages - 1
        nodes_x = x_values[1:step:end, :]
        nodes_s = s_all[1:step:end, :]

        function interpolate(queryPoints)
            out = zeros(length(queryPoints), size(nodes_x, 2))
            for pointIDX = 1:length(queryPoints)
                segment_idx = -1
                for i = 1:length(nodes_s)
                    if nodes_s[i] > queryPoints[pointIDX]
                        segment_idx = i - 1
                        break
                    end
                end

                if segment_idx == -1
                    if queryPoints[pointIDX] < nodes_s[1]
                        segment_idx = 1
                    else
                        segment_idx = length(nodes_s) - 1
                    end
                end
                if segment_idx == 0
                    segment_idx = 1
                end
                segment_start = 1 + (segment_idx - 1) * step
                x_segment = x_values[segment_start:segment_start+step, :]
                u_segment = u_values[segment_start:segment_start+step, :]
                s_segment = s_all[segment_start:segment_start+step]

                derivatives = zeros(stages, size(x_segment, 2))
                for i = 1:stages
                    derivatives[i, :] = f(x_segment[i, :], u_segment[i, :], s_segment[i], false)
                end

                for x_idx = eachindex(x_segment[1, :])
                    polynom = HermiteInterpolation.fit(s_segment, x_segment[:, x_idx], derivatives[:, x_idx])
                    x_interp = polynom(queryPoints[pointIDX])
                    out[pointIDX, x_idx] = x_interp
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

"""
### According to:
    Berrut, J.-P. and Trefethen, L.N. Barycentric Lagrange Interpolation, SIAM Review, 46(3), 2004,
"""
function get_diff_matix(N, variant)
    if variant == "Radau"
        # Left endpoint (τ=-1) + N Radau points (rightmost is τ=+1)
        nodes = reverse([gaussradau(N)[1]..., 1]) .* -1
    elseif variant == "Legendre"
        # Left boundary (τ=-1) + N interior Gauss points + right boundary (τ=+1)
        nodes = [-1.0; gausslegendre(N)[1]; 1.0]
    end

    M = length(nodes)

    # Barycentric weights
    w = ones(M)
    for j in 1:M
        for k in 1:M
            if k != j
                w[j] *= (nodes[j] - nodes[k])
            end
        end
        w[j] = 1.0 / w[j]
    end

    # Full differentiation matrix
    Dfull = zeros(M, M)
    for i in 1:M
        for j in 1:M
            if i != j
                Dfull[i, j] = (w[j] / w[i]) / (nodes[i] - nodes[j])
            end
        end
        Dfull[i, i] = -sum(Dfull[i, :])
    end

    if variant == "Radau"
        # Collocation points are nodes[2:end] (all except left boundary)
        D = Dfull[2:end, :]
    elseif variant == "Legendre"
        # Collocation points are interior nodes only: nodes[2:end-1]
        D = Dfull[2:end-1, :]
    end

    return (D, nodes)
end

function create_gauss_pseudospectral_metod(f, pol_order, variant, model, nControls, nStates, track)
    t0 = time()
    (D, nodes) = get_diff_matix(pol_order, variant)
    println("making diff matrix and nodes: $(round(time() - t0, digits=3))s")
    t0 = time()

    # For Radau:   nodes has pol_order+1 entries, D is pol_order × (pol_order+1)
    # For Legendre: nodes has pol_order+2 entries, D is pol_order × (pol_order+2)
    # In both cases D has pol_order rows (one per collocation point)

    # --- Barycentric weights (Berrut & Trefethen, SIAM Review 2004) ---
    function _bary_weights(pts)
        N = length(pts)
        w = ones(N)
        for j in 1:N, k in 1:N
            k != j && (w[j] *= (pts[j] - pts[k]))
        end
        return 1.0 ./ w
    end

    w_state = _bary_weights(nodes)

    # Control polynomial uses only the interior collocation nodes
    if variant == "Radau"
        ctrl_nodes = nodes[2:end]
    elseif variant == "Legendre"
        ctrl_nodes = nodes[2:end-1]
    end
    w_ctrl = _bary_weights(ctrl_nodes)

    """Evaluate barycentric Lagrange polynomial at τ."""
    function _bary_eval(ref_nodes, weights, values, τ)
        for j in eachindex(ref_nodes)
            abs(τ - ref_nodes[j]) < 100 * eps(Float64) && return Float64(values[j])
        end
        num = sum(weights[j] * values[j] / (τ - ref_nodes[j]) for j in eachindex(ref_nodes))
        den = sum(weights[j]             / (τ - ref_nodes[j]) for j in eachindex(ref_nodes))
        return num / den
    end

    function create_dynamic_constraints(segments, pol_order, initialization)

        # For Radau:   each segment contributes pol_order rows, boundaries shared → segments*pol_order + 1 total
        # For Legendre: each segment contributes pol_order+1 rows, boundaries shared → segments*(pol_order+1) + 1 total
        if variant == "Radau"
            stride     = pol_order          # rows per segment (shared boundary counts once)
            n_seg_nodes = pol_order + 1     # nodes per segment including both endpoints
        elseif variant == "Legendre"
            stride     = pol_order + 1
            n_seg_nodes = pol_order + 2
        end

        nrows = segments * stride + 1
        X = Matrix{VariableRef}(undef, nrows, nStates)
        U = Matrix{VariableRef}(undef, nrows, nControls)

        # FX: collocation residuals, pol_order rows per segment
        FX = Matrix{NonlinearExpr}(undef, pol_order, nStates)

        # Gauss-Legendre quadrature weights (only needed for Legendre endpoint constraint)
        if variant == "Legendre"
            _, w_quad = gausslegendre(pol_order)
            wF = Matrix{NonlinearExpr}(undef, pol_order, nStates)
        end

        for i = 1:nrows
            for j = 1:nStates
                X[i, j] = @variable(model)
            end
            for j = 1:nControls
                U[i, j] = @variable(model)
            end
        end

        segment_edges = collect(LinRange(track.sampleDistances[1], track.sampleDistances[end], segments + 1))
        scaled_nodes  = (nodes .+ 1) / 2   # map [-1,1] → [0,1]

        all_nodes = zeros(nrows)

        for i = 1:segments
            start_idx = (i - 1) * stride + 1
            end_idx   = start_idx + stride          # = i*stride + 1
            h = segment_edges[i+1] - segment_edges[i]

            # Physical positions of all nodes in this segment
            for node = 1:n_seg_nodes
                all_nodes[start_idx + node - 1] = segment_edges[i] + scaled_nodes[node] * h
            end
            seg_nodes = all_nodes[start_idx:end_idx]

            # --- Initial guesses ---
            # Left boundary
            init_boundary = initialization.states(seg_nodes[1])
            for k = 1:nStates
                set_start_value(X[start_idx, k], init_boundary[k])
            end
            init_u_boundary = initialization.controls(seg_nodes[1])
            for k = 1:nControls
                set_start_value(U[start_idx, k], init_u_boundary[k])
            end

            # Interior collocation points: indices start_idx+1 … start_idx+pol_order
            # (for Legendre these are nodes 2…pol_order+1 out of pol_order+2)
            # (for Radau   these are nodes 2…pol_order+1 out of pol_order+1, i.e. all except left)
            for j = 1:pol_order
                node_pos = seg_nodes[j + 1]   # seg_nodes[1]=left boundary, so interior starts at [2]
                init_vals = initialization.states(node_pos)
                for k = 1:nStates
                    set_start_value(X[start_idx + j, k], init_vals[k])
                end
                init_u = initialization.controls(node_pos)
                for k = 1:nControls
                    set_start_value(U[start_idx + j, k], init_u[k])
                end
                # Evaluate dynamics at collocation point
                FX[j, :] = f(X[start_idx + j, :], U[start_idx + j, :], node_pos)
            end

            # For Legendre: also initialise the right boundary node
            if variant == "Legendre"
                right_pos = seg_nodes[end]   # = seg_nodes[pol_order+2]
                init_right = initialization.states(right_pos)
                for k = 1:nStates
                    set_start_value(X[end_idx, k], init_right[k])
                end
                init_u_right = initialization.controls(right_pos)
                for k = 1:nControls
                    set_start_value(U[end_idx, k], init_u_right[k])
                end
            end

            # --- Differentiation matrix constraint (enforces ODE at collocation points) ---
            # D * X_seg = (h/2) * F  where X_seg has n_seg_nodes rows
            @constraint(model, D * X[start_idx:end_idx, :] .== (h / 2) .* FX)

            # --- Legendre only: quadrature constraint to pin the right endpoint ---
            # X(τ=+1) = X(τ=-1) + (h/2) * Σ w_i f(x_i, u_i)
            if variant == "Legendre"
                for j = 1:pol_order
                    wF[j, :] = FX[j, :]   # reuse already-evaluated FX
                end
                @constraint(model, X[end_idx, :] .== X[start_idx, :] + (h / 2) * wF' * w_quad)
            end
        end

        println("making constraints: $(round(time() - t0, digits=3))s")
        return [model, X, U, all_nodes, segment_edges]
    end

    """
    Build interpolating functions from the optimised solution.

    Per Garg et al. (2010):
    - States:   degree-N (Radau) or degree-(N+1) (Legendre) Lagrange polynomial through all nodes.
    - Controls: degree-(N-1) polynomial through interior collocation nodes only.
    Barycentric form per Berrut & Trefethen (2004).
    """
    function create_gauss_interpolator(X_vals, U_vals, all_nodes, segment_edges)
        segments = length(segment_edges) - 1

        if variant == "Radau"
            stride = pol_order
        elseif variant == "Legendre"
            stride = pol_order + 1
        end

        function find_segment(s)
            for i in 1:segments
                s <= segment_edges[i+1] && return i
            end
            return segments
        end

        function state_interp(s)
            seg = find_segment(s)
            h   = segment_edges[seg+1] - segment_edges[seg]
            τ   = 2 * (s - segment_edges[seg]) / h - 1
            i0  = (seg - 1) * stride + 1
            x_seg = X_vals[i0 : i0 + stride, :]          # (n_seg_nodes) × nStates
            return [_bary_eval(nodes, w_state, x_seg[:, k], τ) for k in 1:size(x_seg, 2)]
        end

        function control_interp(s)
            seg = find_segment(s)
            h   = segment_edges[seg+1] - segment_edges[seg]
            τ   = 2 * (s - segment_edges[seg]) / h - 1
            i0  = (seg - 1) * stride + 1
            # Interior collocation nodes only (exclude boundary/ies)
            u_seg = U_vals[i0+1 : i0 + pol_order, :]     # pol_order × nControls

            hold_controls = 1
            if hold_controls == 1
                u0 = u_seg[1, :]
                return [u0[k] for k in 1:size(u_seg, 2)]
            else
                return [_bary_eval(ctrl_nodes, w_ctrl, u_seg[:, k], τ) for k in 1:size(u_seg, 2)]
            end
        end

        return Result_interpolation(state_interp, control_interp, all_nodes)
    end

    GaussMethod = Collocation(
        create_dynamic_constraints,
        create_gauss_interpolator,
        nothing,
        nothing,
        f
    )
    return GaussMethod
end