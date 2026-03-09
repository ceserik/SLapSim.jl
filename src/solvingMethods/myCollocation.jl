using Revise
using RungeKutta
using JuMP, Ipopt, Zygote
using Infiltrator
using HermiteInterpolation
using FastGaussQuadrature
using Symbolics
using Nemo



struct Collocation
    createConstraints
    createInterpolator
    tableau
    interpolator::Union{Function,Nothing}
    f
end

# This function creates Legandre-Gauss-Lobatto Collocation
function createLobattoIIIA(stage, f)
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

function Lagrange_polynomial_basis(nodes::Vector{Float64},τ::Num)
    k = length(nodes)
    #println(nodes)
    
    L = zero(Num)
    l = ones(Num, k)
    for j = 1:k
        l[j] = one(Num)
        for m = 1:k
            if j !=m
                l[j] *= (τ - nodes[m])/(nodes[j] - nodes[m])
            end
        end
        #L += nodes[j]*l[j]
    end
    return l
end

function get_diff_matix(N,variant)
    if variant == "Radau"
        nodes = reverse([gaussradau(N)[1]..., 1]) .* -1
        D = zeros(N, N + 1)
    end

    if variant == "Legendre"
        nodes = [-1; gausslegendre(N)[1]]
        D = zeros(N, N + 1)
    end

    @Symbolics.variables τ
    L = Lagrange_polynomial_basis(nodes, τ)
    L_dot = expand(Symbolics.derivative(L, τ))

    for i = 1:N
        D_1 = substitute(L_dot, Dict(τ => nodes[i+1]))
        D[i, :] = Symbolics.symbolic_to_float.(D_1)
    end
    return (D, nodes)
end

function create_gauss_pseudospectral_metod(f,pol_order,variant,model,nControls,nStates,track)
    
    (D, nodes) = get_diff_matix(pol_order, variant)
    FX = Matrix{NonlinearExpr}(undef,  pol_order , nStates)
    function create_dynamic_constraints(segments,pol_order,initialization)

        X = Matrix{VariableRef}(undef, segments * pol_order + 1, nStates)
        U = Matrix{VariableRef}(undef, segments * pol_order + 1, nControls)

        for i = 1:segments * pol_order +1
            for j = 1:nStates
                X[i, j] = @variable(model)    
            end
            
            
            for j = 1:nControls
                U[i, j] = @variable(model)
            end
        end
        
        segment_edges = LinRange(track.sampleDistances[1],track.sampleDistances[end],segments+1)
        scaled_nodes = (nodes.+1)/2
        segment_nodes = zeros(length(nodes))
        for i = 1:length(segment_edges)-1
            start_idx = (i-1)*pol_order+1
            end_idx = i*pol_order +1
            h = segment_edges[i+1] - segment_edges[i]
            
            for node = 1:length(nodes)
                segment_nodes[node] = segment_edges[i] + scaled_nodes[node]*(segment_edges[i+1] - segment_edges[i])
            end

            # Initialize boundary point
            init_boundary = initialization.states(segment_nodes[1])
            for k = 1:nStates
                set_start_value(X[start_idx, k], init_boundary[k])
            end
            init_u_boundary = initialization.controls(segment_nodes[1])
            for k = 1:nControls
                set_start_value(U[start_idx, k], init_u_boundary[k])
            end

            for j = 1:pol_order
                # j-th collocation point corresponds to nodes[j+1]
                init_vals = initialization.states(segment_nodes[j+1])
                for k = 1:nStates
                    set_start_value(X[start_idx+j, k], init_vals[k])
                end
                init_u = initialization.controls(segment_nodes[j+1])
                for k = 1:nControls
                    set_start_value(U[start_idx+j, k], init_u[k])
                end
                FX[j,:] = f(X[start_idx+j,:], U[start_idx+j,:], segment_nodes[j+1])
            end
            @constraint(model, D*X[start_idx:end_idx,:] .== (h/2) .* FX)
        end

        return [model, X, U, segment_nodes]
    end
    GaussMethod = Collocation(
        create_dynamic_constraints,
        nothing,
        nothing,
        nothing,
        f
    )
    return GaussMethod
end


