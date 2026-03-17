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



"""
Darby, Christopher L., William W. Hager, and Anil V. Rao. 
“An Hp ‐adaptive Pseudospectral Method for Solving Optimal Control Problems.”
Optimal Control Applications and Methods 32, no. 4 (2011): 476–502.
https://doi.org/10.1002/oca.957.

"""
function create_gauss_legendre(f, pol_order, variant, model, nControls, nStates, track)
    (nodes_LG,w2) = gausslegendre(pol_order)
    τ = [-1; nodes_LG]
    w = barycentric_weights(τ)

    D = diff_matrix(τ, w)
    D = D[2:end, :]

    function create_dynamic_constraints(segments, initialization)
        X = Matrix{VariableRef}(undef, segments * (pol_order + 1)+1, nStates)
        U = Matrix{VariableRef}(undef, segments * (pol_order), nControls)
        FX = Matrix{NonlinearExpr}(undef, pol_order, nStates)

        #this creation of cariables is very clumsy, buut when they are crated interlaced(like this), the hessian an lagrangian are diagonal banded matrices
        for i = 1:segments*(pol_order+1)+1
            for j = 1:nStates
                X[i, j] = @variable(model)
            end

            if mod(i - 1, pol_order + 1) != 0
                for j = 1:nControls
                    seg = div(i - 1, pol_order + 1)
                    k = mod(i - 1, pol_order + 1)
                    U[seg*pol_order+k, j] = @variable(model)
                end
            end
        end

        segment_edges = LinRange(track.sampleDistances[1], track.sampleDistances[end], segments+1)
        x_start_idx = 1
        u_start_idx = 1
        scaled_τ = (τ .+ 1) / 2 #difference between nodes_LG and tau is that tau has 1 extra point which is not collocated -1
        all_s = zeros(segments * (pol_order + 1))

        #now create the collocation constraints
        for i = 1:segments
            end_idx = x_start_idx + pol_order
            
            h = segment_edges[i+1] - segment_edges[i]

            for node = 1:length(τ)
                all_s[x_start_idx+node-1] = segment_edges[i] + scaled_τ[node] * h
            end
            seg_nodes = all_s[x_start_idx:end_idx]
            
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
                
            end
#            @infiltrate
            @constraint(model, D * X[x_start_idx:end_idx, :] .== (h / 2) .* FX)
            @constraint(model,X[end_idx+1,:] .== X[x_start_idx,:] + (h / 2) * FX' * w2)
            x_start_idx  = end_idx + 1
            #u_start_idx += pol_order 

            
        end
        return [model, X, U, all_s, segment_edges]
    end


    GaussMethod = Collocation(
        create_dynamic_constraints,
        nothing,
        nothing,
        nothing,
        nothing
    )
end




#track = doubleTurn(false, 0.5)
#car = createSimplestSingleTrack()
#model = JuMP.Model(Ipopt.Optimizer)
#initialization = initializeSolution_interpolation(car, track, 200)
#xd = create_gauss_legendre(1, 2, 2, model, 1, 1, track)
#xd.createConstraints(2, initialization)