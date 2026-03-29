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
    #tableau = TableauRunge()
    stages = tableau.s
    function createDynamicConstraints(f, Xsize::Int64, Usize::Int64, iterpolationFunction, sampplingPoints::AbstractVector{Float64}, model::JuMP.Model, X_init, U_init)
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

        # Bounds: X = [vx, vy, ψ, ψ̇, n, t], U = [torque, steering]
        

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
                h_node = node_s[i+1] - node_s[i]
                for st = 2:(stages-1)
                    s_all[idx] = node_s[i] + h_node * tableau.c[st]
                    idx += 1
                end
            end
        end
        constrained = falses(totalPoints)
        # Collocation / continuity constraints for each interval
        for interval = 1:(N-1)
            interval_start_idx = 1 + (interval - 1) * (stages - 1)
            node_start_s = node_s[interval]
            node_end_s = node_s[interval+1]
            h = node_end_s - node_start_s

            for stage = 2:stages
                # index in X/U corresponding to this collocation/stage point
                next_idx = interval_start_idx + stage - 1


                # Build Runge-Kutta weighted sum of f evaluated at the stage points in this interval
                fxSum = zeros(NonlinearExpr, Xsize)
                for col = 1:stages
                    x_idx = interval_start_idx + col - 1
                    fxSum += tableau.a[stage, col] * f(X[x_idx, :], U[x_idx, :], s_all[x_idx], nothing)
                    #xdddd = f(X[x_idx, :], U[x_idx, :], s_all[x_idx],true)
                end

                # Only add model constraints for each time point once
                if !constrained[next_idx]
                    _ = f(X[next_idx, :], U[next_idx, :], s_all[next_idx], model)
                    constrained[next_idx] = true
                end
                if !constrained[interval_start_idx]
                    _ = f(X[interval_start_idx, :], U[interval_start_idx, :], s_all[interval_start_idx], model)
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


    function createLobattoInterpolator(x_values, u_values, s_all)
        step = stages - 1
        nodes_x = x_values[1:step:end, :]
        nodes_s = s_all[1:step:end, :]

        function interpolate(queryPoints)
            out = zeros(length(queryPoints), size(nodes_x, 2))
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
        nodes = reverse([gaussradau(N)[1]..., 1]) .* -1
    elseif variant == "Legendre"
        nodes = [-1; gausslegendre(N)[1]]
    elseif variant == "Lobatto"
        nodes = gausslobatto(N)[1]
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

    # Extract rows for collocation points only (nodes[2:end])
    # Skip row 1 because the left boundary is constrained by continuity, not as an ODE point
    D = Dfull[2:end, :]
    return (D, nodes)
end

function create_gauss_pseudospectral_metod(f, pol_order, variant, model, nControls, nStates, track)
    t0 = time()
    (D, nodes) = get_diff_matix(pol_order, variant)
    println("making diff matrix and nodes: $(round(time() - t0, digits=3))s")
    t0 = time()


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
    
    ctrl_nodes = nodes[2:end]
    w_ctrl = _bary_weights(ctrl_nodes)

    """Evaluate barycentric Lagrange polynomial at τ given reference nodes, weights, values."""
    function _bary_eval(ref_nodes, weights, values, τ)
        for j in eachindex(ref_nodes)
            abs(τ - ref_nodes[j]) < 100 * eps(Float64) && return Float64(values[j])
        end
        num = sum(weights[j] * values[j] / (τ - ref_nodes[j]) for j in eachindex(ref_nodes))
        den = sum(weights[j] / (τ - ref_nodes[j]) for j in eachindex(ref_nodes))
        return num / den
    end

    function create_dynamic_constraints(segments, initialization)

        if variant == "Radau"
            additional_nodes = 0
            final = 1
        elseif variant == "Legendre"
            additional_nodes = 1
            final = 1
        elseif variant == "Lobatto"
            additional_nodes = -1
            final = 0
        end

        # X allocation: unchanged from original
        if variant == "Radau" || variant == "Legendre"
            X = Matrix{VariableRef}(undef, segments * (pol_order + additional_nodes) + final, nStates)
        else
            X = Matrix{VariableRef}(undef, (segments - 1) * (pol_order - 1) + pol_order, nStates)
        end


        if variant == "Legendre"
            U = Matrix{VariableRef}(undef, segments * pol_order, nControls)
        elseif variant == "Radau"
            U = Matrix{VariableRef}(undef, segments * (pol_order + additional_nodes) + final, nControls)
        else
            U = Matrix{VariableRef}(undef, (segments - 1) * (pol_order - 1) + pol_order, nControls)
        end

        n_ode_rows = size(D, 1)
        FX = Matrix{NonlinearExpr}(undef, n_ode_rows, nStates)
        if variant == "Legendre"
            kl, w = gausslegendre(pol_order)
        end

        # Bounds: X = [vx, vy, ψ, ψ̇, n, t], U = [torque, steering]
        x_lb = [0.5, -20.0, -2π, -5.0, -10.0, 0.0]
        x_ub = [40.0, 20.0,  2π,  5.0,  10.0, 200.0]
        u_lb = [-29.0, -20/180*π, -29.0]
        u_ub = [ 29.0,  20/180*π,  29.0]

        # Create X variables
        for i = 1:size(X, 1)
            for j = 1:nStates
                X[i, j] = @variable(model, lower_bound=x_lb[j], upper_bound=x_ub[j])
            end
        end

        for i = 1:size(U, 1)
            for j = 1:nControls
                U[i, j] = @variable(model, lower_bound=u_lb[j], upper_bound=u_ub[j])
            end
        end

        segment_edges = collect(LinRange(track.sampleDistances[1], track.sampleDistances[end], segments + 1))
        scaled_nodes = (nodes .+ 1) / 2

        all_nodes = zeros(segments * (pol_order + additional_nodes) + 1)

        x_start_idx = 1
        u_start_idx = 1 

        for i = 1:segments
            if variant == "Legendre"
                end_idx = x_start_idx + pol_order
            elseif variant == "Radau"
                end_idx = x_start_idx + pol_order
            elseif variant == "Lobatto"
                end_idx = x_start_idx + pol_order - 1
            end

            h = segment_edges[i+1] - segment_edges[i]

            for node = 1:length(nodes)
                all_nodes[x_start_idx + node - 1] = segment_edges[i] + scaled_nodes[node] * h
            end
            seg_nodes = all_nodes[x_start_idx:end_idx]

            
            init_boundary = initialization.states(seg_nodes[1])
            for k = 1:nStates
                set_start_value(X[x_start_idx, k], init_boundary[k])
            end

            
            if variant != "Legendre"
                init_u_boundary = initialization.controls(seg_nodes[1])
                for k = 1:nControls
                    set_start_value(U[x_start_idx, k], init_u_boundary[k])
                end
            end

            for j = 1:n_ode_rows
                node_s = seg_nodes[j + 1]   
                x_idx  = x_start_idx + j
                
                u_idx  = (variant == "Legendre") ? u_start_idx + j - 1 : x_idx

                init_vals = initialization.states(node_s)
                for k = 1:nStates
                    set_start_value(X[x_idx, k], init_vals[k])
                end
                init_u = initialization.controls(node_s)
                for k = 1:nControls
                    set_start_value(U[u_idx, k], init_u[k])
                end

                FX[j, :] = f(X[x_idx, :], U[u_idx, :], node_s)  # FIX: use u_idx
            end

            @constraint(model, D * X[x_start_idx:end_idx, :] .== (h / 2) .* FX)

            if variant == "Legendre"
                @constraint(model, X[end_idx+1, :] .== X[x_start_idx, :] + (h / 2) * FX' * w)
                all_nodes[end_idx+1] = segment_edges[i+1]
                x_start_idx  = end_idx + 1
                u_start_idx += pol_order   # FIX: advance by pol_order (no boundary slot)
            elseif variant == "Radau"
                x_start_idx = end_idx
            else
                x_start_idx = end_idx
            end
        end

        println("making constraints: $(round(time() - t0, digits=3))s")

        return [model, X, U, all_nodes, segment_edges]
    end

    """
    Build interpolating functions from the optimised solution using the LGR polynomial.

    Per Garg et al. (2010) "A Unified Framework for the Numerical Solution of Optimal Control
    Problems Using Pseudospectral Methods":
    - States are represented by a degree-N Lagrange polynomial through all N+1 nodes
      {τ₀=-1, τ₁,…,τ_N} (initial boundary + N LGR collocation points).
    - Controls are represented by a degree-(N-1) Lagrange polynomial through the N LGR
      collocation nodes {τ₁,…,τ_N} only.
    Evaluation uses the numerically stable barycentric form (Berrut & Trefethen 2004).
    """
    function create_gauss_interpolator(X_vals, U_vals, all_nodes, segment_edges; hold_controls=true)
        segments = length(segment_edges) - 1

        x_stride = variant == "Legendre" ? pol_order + 1 :
                   variant == "Lobatto"  ? pol_order - 1 : pol_order

        u_stride = variant == "Legendre" ? pol_order : x_stride

        function find_segment(s)
            for i in 1:segments
                s <= segment_edges[i+1] && return i
            end
            return segments
        end

        function state_interp(s)
            seg   = find_segment(s)
            h     = segment_edges[seg+1] - segment_edges[seg]
            τ     = 2 * (s - segment_edges[seg]) / h - 1
            i0    = (seg - 1) * x_stride + 1
            i_end = min(i0 + pol_order, size(X_vals, 1))
            x_seg = X_vals[i0:i_end, :]
            return [_bary_eval(nodes, w_state, x_seg[:, k], τ) for k in 1:size(x_seg, 2)]
        end

        u_positions = if variant == "Legendre"
            pos = zeros(segments * pol_order)
            for i in 1:segments
                x_start = (i - 1) * x_stride + 1
                u_base  = (i - 1) * pol_order
                for j in 1:pol_order
                    pos[u_base + j] = all_nodes[x_start + j]
                end
            end
            pos
        else
            all_nodes
        end

        function control_interp(s)
            seg   = find_segment(s)
            h     = segment_edges[seg+1] - segment_edges[seg]
            τ     = 2 * (s - segment_edges[seg]) / h - 1
            i0    = (seg - 1) * u_stride + 1
            i_end = i0 + u_stride - 1
            u_seg = U_vals[i0:i_end, :]
            return [_bary_eval(ctrl_nodes, w_ctrl, u_seg[:, k], τ) for k in 1:size(u_seg, 2)]
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