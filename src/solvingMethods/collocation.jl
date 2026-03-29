using Revise
using JuMP
using Infiltrator
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
function create_gauss_legendre(f, pol_order, variant, model, nControls, nStates, track)

    (nodes_LG, w2) = gausslegendre(pol_order)
    τ = [-1; nodes_LG]
    w = barycentric_weights(τ)

    D = diff_matrix(τ, w)
    D = D[2:end, :]

    function create_dynamic_constraints(segments, initialization)
        X = Matrix{VariableRef}(undef, segments * (pol_order + 1) + 1, nStates)
        U = Matrix{VariableRef}(undef, segments * (pol_order), nControls)
        FX = Matrix{NonlinearExpr}(undef, pol_order, nStates)

        #this creation of cariables is very clumsy, buut when they are crated interlaced(like this), the hessian an lagrangian are diagonal banded matrices
        # Bounds: X = [vx, vy, ψ, ψ̇, n, t], U = [torque, steering]
        x_lb = [0.5, -20.0, -2π, -5.0, -10.0, 0.0]
        x_ub = [40.0, 20.0,  2π,  5.0,  10.0, 200.0]
        u_lb = [-29.0, -20/180*π, -29.0]
        u_ub = [ 29.0,  20/180*π,  29.0]

        for i = 1:segments*(pol_order+1)+1
            for j = 1:nStates
                X[i, j] = @variable(model, lower_bound=x_lb[j], upper_bound=x_ub[j])
            end

            if mod(i - 1, pol_order + 1) != 0
                for j = 1:nControls
                    seg = div(i - 1, pol_order + 1)
                    k = mod(i - 1, pol_order + 1)
                    U[seg*pol_order+k, j] = @variable(model, lower_bound=u_lb[j], upper_bound=u_ub[j])
                end
            end
        end

        segment_edges = LinRange(track.sampleDistances[1], track.sampleDistances[end], segments + 1)
        x_start_idx = 1
        u_start_idx = 1
        scaled_τ = (τ .+ 1) / 2 #difference between nodes_LG and tau is that tau has 1 extra point which is not collocated -1
        all_s = zeros(segments * (pol_order + 1))

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
                u_start_idx += 1

                init_vals = initialization.states(node_s)
                for k = 1:nStates
                    set_start_value(X[x_idx, k], init_vals[k])
                end

                init_u = initialization.controls(node_s)
                for k = 1:nControls
                    set_start_value(U[u_idx, k], init_u[k])
                end

                FX[j, :] = f(X[x_idx, :], U[u_idx, :], node_s)
                last_u_idx = u_idx
            end
            init_boundary = initialization.states(segment_edges[i+1])
            for k = 1:nStates
                set_start_value(X[end_idx+1, k], init_boundary[k])
            end
            @constraint(model, D * X[x_start_idx:end_idx, :] .== (h / 2) .* FX)
            @constraint(model, X[end_idx+1, :] .== X[x_start_idx, :] + (h / 2) * FX' * w2)
            #last_u_idx = u_start_idx - 1
            f(X[end_idx+1, :], U[last_u_idx, :], segment_edges[i+1])
            
            x_start_idx = end_idx + 1
            #u_start_idx += pol_order 


        end
        println("Creating constraints: $(round(time() - t0_constraints, digits=3))s")
        return [model, X, U, all_s, segment_edges]
    end

    function create_interpolation(X_vals, U_vals, all_s, segment_edges)


        x_stride = pol_order + 1
        u_stride = pol_order 
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

        return Result_interpolation(state_interp, control_interp, all_s)
    end

    GaussMethod = Collocation(
        create_dynamic_constraints,
        create_interpolation,
        nothing,
        nothing,
        nothing
    )
end

