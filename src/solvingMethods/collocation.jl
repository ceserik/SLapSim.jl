using JuMP

using FastGaussQuadrature
using SLapSim
#include("../solvingMethods/optInterface.jl")
#struct Collocation
#    createConstraints
#    createInterpolator
#    tableau
#    interpolator::Union{Function,Nothing}
#    f
#end

"""
Berrut, Jean-Paul, and Lloyd N. Trefethen. 
“Barycentric Lagrange Interpolation.” 
SIAM Review 46, no. 3 (2004): 501–17.
https://doi.org/10.1137/S0036144502417715.
"""
function barycentric_weights(x)
    w = zeros(length(x))

    for j = eachindex(x)
        tmp = 1
        for k = eachindex(x)
            if k != j
                tmp *= (x[j] - x[k])
            end
        end
        w[j] = 1 / tmp
    end
    return w
end


function diff_matrix(x, w)
    D = zeros(length(x), length(x))
    for i = eachindex(x)
        for j = eachindex(x)
            if j != i
                D[i, j] = (w[j] / w[i]) / (x[i] - x[j])
            end
        end
    end

    for i = eachindex(w)
        sum = 0
        for j = eachindex(w)
            if j != i
                sum -= D[i, j]
            end
        end
        D[i, i] = sum

    end
    return D
end


function _bary_eval(ref_nodes, weights, values, τ)
    for j in eachindex(ref_nodes)
        abs(τ - ref_nodes[j]) < 100 * eps(Float64) && return Float64(values[j])
    end
    num = sum(weights[j] * values[j] / (τ - ref_nodes[j]) for j in eachindex(ref_nodes))
    den = sum(weights[j] / (τ - ref_nodes[j]) for j in eachindex(ref_nodes))
    return num / den
end


"""
Darby, Christopher L., William W. Hager, and Anil V. Rao. 
“An Hp ‐adaptive Pseudospectral Method for Solving Optimal Control Problems.”
Optimal Control Applications and Methods 32, no. 4 (2011): 476–502.
https://doi.org/10.1002/oca.957.

"""
function create_gauss_legendre(f, pol_order, variant, model, nControls, nStates, track; car=nothing)

    (nodes_LG, w2) = gausslegendre(pol_order)
    τ = [-1; nodes_LG]
    w = barycentric_weights(τ)

    D = diff_matrix(τ, w)
    D = D[2:end, :]

    function create_dynamic_constraints(segment_edges, initialization)
        segments = length(segment_edges) - 1
        X = Matrix{VariableRef}(undef, segments * (pol_order + 1) + 1, nStates)
        U = Matrix{VariableRef}(undef, segments * (pol_order), nControls)
        FX = Matrix{NonlinearExpr}(undef, pol_order, nStates)

        # Placeholder scaled bounds; real bounds applied later by applyBounds!.
        for i = 1:segments*(pol_order+1)+1
            for j = 1:nStates
                X[i, j] = @variable(model, lower_bound=-1.0, upper_bound=1.0)
            end

            if mod(i - 1, pol_order + 1) != 0
                for j = 1:nControls
                    seg = div(i - 1, pol_order + 1)
                    k = mod(i - 1, pol_order + 1)
                    U[seg*pol_order+k, j] = @variable(model, lower_bound=-1.0, upper_bound=1.0)
                end
            end
        end

        x_start_idx = 1
        u_start_idx = 1
        scaled_τ = (τ .+ 1) / 2 #difference between nodes_LG and tau is that tau has 1 extra point which is not collocated -1
        all_s = zeros(segments * (pol_order + 1) + 1)
        s_controls = zeros(segments * pol_order)

        #now create the collocation constraints
        t0_constraints = time()
        for i = 1:segments
            end_idx = x_start_idx + pol_order

            h = segment_edges[i+1] - segment_edges[i]

            for node = 1:length(τ)
                all_s[x_start_idx+node-1] = segment_edges[i] + scaled_τ[node] * h
            end
            seg_nodes = all_s[x_start_idx:end_idx]
            last_u_idx = -1
            for j = 1:length(nodes_LG)
                node_s = seg_nodes[j+1]
                x_idx = x_start_idx + j
                u_idx = u_start_idx
                s_controls[u_idx] = node_s
                u_start_idx += 1

                init_vals = initialization.states(node_s)
                for k = 1:nStates
                    set_start_value(X[x_idx, k], init_vals[k])
                end

                init_u = initialization.controls(node_s)
                for k = 1:nControls
                    set_start_value(U[u_idx, k], init_u[k])
                end

                FX[j, :] = f(X[x_idx, :], U[u_idx, :], node_s, model)
                last_u_idx = u_idx
            end
            init_boundary = initialization.states(segment_edges[i+1])
            for k = 1:nStates
                set_start_value(X[end_idx+1, k], init_boundary[k])
            end
            @constraint(model, D * X[x_start_idx:end_idx, :] .== (h / 2) .* FX)
            @constraint(model, X[end_idx+1, :] .== X[x_start_idx, :] + (h / 2) * FX' * w2)
            #last_u_idx = u_start_idx - 1
            f(X[end_idx+1, :], U[last_u_idx, :], segment_edges[i+1], model)
            
            x_start_idx = end_idx + 1
            #u_start_idx += pol_order 


        end
        all_s[end] = segment_edges[end]
        println("Creating constraints: $(round(time() - t0_constraints, digits=3))s")
        return [model, X, U, all_s, segment_edges, s_controls]
    end

    function create_interpolation(X_vals, U_vals, all_s, segment_edges)
        segments = length(segment_edges) - 1

        x_stride = pol_order + 1
        u_stride = pol_order
        function find_segment(s)
            for i in 1:segments
                s < segment_edges[i+1] && return i
            end
            return segments
        end

        function state_interp(s)
            seg = find_segment(s)
            h = segment_edges[seg+1] - segment_edges[seg]
            τ_eval = 2 * (s - segment_edges[seg]) / h - 1   
            i0 = (seg - 1) * x_stride + 1
            i_end = min(i0 + pol_order, size(X_vals, 1))
            x_seg = X_vals[i0:i_end, :]
            return [_bary_eval(τ, w, x_seg[:, k], τ_eval) for k in 1:size(x_seg, 2)]
        end

        w_u = barycentric_weights(nodes_LG)

        function control_interp(s)
            seg = find_segment(s)
            h = segment_edges[seg+1] - segment_edges[seg]
            τ_eval = 2 * (s - segment_edges[seg]) / h - 1
            i0 = (seg - 1) * u_stride + 1
            i_end = i0 + pol_order - 1
            u_seg = U_vals[i0:i_end, :]
            return [_bary_eval(nodes_LG, w_u, u_seg[:, k], τ_eval) for k in 1:size(u_seg, 2)]
        end

        return Result_interpolation(state_interp, control_interp, all_s, collect(Float64, segment_edges))
    end

    GaussMethod = Collocation(
        create_dynamic_constraints,
        create_interpolation,
        nothing,
        nothing,
        f
    )
    return GaussMethod
end


"""
Garg, D. et al. "A Unified Framework for the Numerical Solution of Optimal Control Problems
Using Pseudospectral Methods." Automatica 46 (2010): 1843–1851.

Flipped Legendre-Gauss-Radau (LGR): collocation includes the right endpoint (+1).
nodes layout per segment: [-1 (non-collocated left), τ₁, …, τ_{N-1}, +1 (collocated right)].
Adjacent segments share the right/left endpoint variable, giving free state continuity.
"""
function create_radau(f, pol_order, model, nControls, nStates, track)
    # Flipped LGR: gaussradau gives N nodes incl. -1; append +1, reverse, negate.
    nodes = reverse([gaussradau(pol_order)[1]..., 1]) .* -1
    w = barycentric_weights(nodes)
    Dfull = diff_matrix(nodes, w)
    D = Dfull[2:end, :]  # collocation rows only

    function create_dynamic_constraints(segment_edges, initialization)
        segments = length(segment_edges) - 1
        nPoints = segments * pol_order + 1
        X = Matrix{VariableRef}(undef, nPoints, nStates)
        U = Matrix{VariableRef}(undef, nPoints, nControls)
        FX = Matrix{NonlinearExpr}(undef, pol_order, nStates)
        all_s = zeros(nPoints)

        # Placeholder scaled bounds; real bounds applied later by applyBounds!.
        for i = 1:nPoints
            for j = 1:nStates
                X[i, j] = @variable(model, lower_bound=-1.0, upper_bound=1.0)
            end
            for j = 1:nControls
                U[i, j] = @variable(model, lower_bound=-1.0, upper_bound=1.0)
            end
        end

        x_start_idx = 1
        scaled_τ = (nodes .+ 1) / 2

        t0_constraints = time()
        for i = 1:segments
            end_idx = x_start_idx + pol_order
            h = segment_edges[i+1] - segment_edges[i]

            for node = 1:length(nodes)
                all_s[x_start_idx + node - 1] = segment_edges[i] + scaled_τ[node] * h
            end
            seg_nodes = all_s[x_start_idx:end_idx]

            # Init non-collocated left boundary (state + control).
            init_left   = initialization.states(seg_nodes[1])
            init_u_left = initialization.controls(seg_nodes[1])
            for k = 1:nStates;   set_start_value(X[x_start_idx, k], init_left[k]);   end
            for k = 1:nControls; set_start_value(U[x_start_idx, k], init_u_left[k]); end

            # Collocation points: nodes[2..end] map to indices x_start_idx+1..end_idx.
            for j = 1:pol_order
                node_s = seg_nodes[j+1]
                idx = x_start_idx + j

                init_vals = initialization.states(node_s)
                for k = 1:nStates;   set_start_value(X[idx, k], init_vals[k]); end
                init_u = initialization.controls(node_s)
                for k = 1:nControls; set_start_value(U[idx, k], init_u[k]);    end

                FX[j, :] = f(X[idx, :], U[idx, :], node_s, model)
            end

            @constraint(model, D * X[x_start_idx:end_idx, :] .== (h / 2) .* FX)

            # Adjacent segments share the right endpoint (= next segment's left bdry).
            x_start_idx = end_idx
        end
        # Pin the very first (non-collocated) control to the first collocation point
        # so it isn't a free variable causing garbage at s=segment_edges[1].
        @constraint(model, U[1, :] .== U[2, :])

        # C¹ and C² control continuity across segment boundaries:
        # k-th derivative from left segment at τ=+1 == k-th derivative from right segment at τ=-1
        # (chain rule factor (2/h)^k folds the [-1,1]→[a,b] mapping).
        #D2 = Dfull * Dfull
        #d_right  = Dfull[end, :]   # ∂/∂τ at τ=+1
        #d_left   = Dfull[1, :]     # ∂/∂τ at τ=-1
        #d2_right = D2[end, :]      # ∂²/∂τ² at τ=+1
        #d2_left  = D2[1, :]        # ∂²/∂τ² at τ=-1
        #for seg = 1:segments-1
        #    i_left  = (seg - 1) * pol_order + 1
        #    i_right = seg * pol_order + 1
        #    h_l = segment_edges[seg+1] - segment_edges[seg]
        #    h_r = segment_edges[seg+2] - segment_edges[seg+1]
        #    for k = 1:nControls
        #        @constraint(model,
        #            (2/h_l) * sum(d_right[j] * U[i_left  + j - 1, k] for j = 1:pol_order+1) ==
        #            (2/h_r) * sum(d_left[j]  * U[i_right + j - 1, k] for j = 1:pol_order+1))
        #        @constraint(model,
        #            (2/h_l)^2 * sum(d2_right[j] * U[i_left  + j - 1, k] for j = 1:pol_order+1) ==
        #            (2/h_r)^2 * sum(d2_left[j]  * U[i_right + j - 1, k] for j = 1:pol_order+1))
        #    end
        #end

        println("Creating constraints: $(round(time() - t0_constraints, digits=3))s")
        # U lives on the same nodes as X → s_controls == all_s.
        return [model, X, U, all_s, segment_edges, all_s]
    end

    function create_interpolation(X_vals, U_vals, all_s, segment_edges)
        segments = length(segment_edges) - 1

        function find_segment(s)
            for i in 1:segments
                s <= segment_edges[i+1] && return i
            end
            return segments
        end

        function state_interp(s)
            seg = find_segment(s)
            h = segment_edges[seg+1] - segment_edges[seg]
            τ_eval = 2 * (s - segment_edges[seg]) / h - 1
            i0 = (seg - 1) * pol_order + 1
            x_seg = X_vals[i0 : i0 + pol_order, :]
            return [_bary_eval(nodes, w, x_seg[:, k], τ_eval) for k in 1:size(x_seg, 2)]
        end

        function control_interp(s)
            seg = find_segment(s)
            h = segment_edges[seg+1] - segment_edges[seg]
            τ_eval = 2 * (s - segment_edges[seg]) / h - 1
            i0 = (seg - 1) * pol_order + 1
            u_seg = U_vals[i0 : i0 + pol_order, :]
            return [_bary_eval(nodes, w, u_seg[:, k], τ_eval) for k in 1:size(u_seg, 2)]
        end

        return Result_interpolation(state_interp, control_interp, all_s, collect(Float64, segment_edges))
    end

    return Collocation(
        create_dynamic_constraints,
        create_interpolation,
        nothing,
        nothing,
        f
    )
end

