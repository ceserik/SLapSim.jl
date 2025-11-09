using JuMP
import Ipopt
import LinearAlgebra
import Random
import SparseArrays
import MathOptInterface as MOI
using UnicodePlots

function fill_off_diagonal(H)
    ret = H + H'
    row_vals = SparseArrays.rowvals(ret)
    non_zeros = SparseArrays.nonzeros(ret)
    for col in 1:size(ret, 2)
        for i in SparseArrays.nzrange(ret, col)
            if col == row_vals[i]
                non_zeros[i] /= 2
            end
        end
    end
    return ret
end





function compute_optimal_hessian(model::Model)
    rows = Any[]
    nlp = MOI.Nonlinear.Model()
    
    # Get constraints directly from the backend in index order
    for (F, S) in list_of_constraint_types(model)
        if F == VariableRef || S <: MOI.AbstractScalarSet && F == VariableRef
            continue  # Skip variable bounds
        end
        cons = all_constraints(model, F, S)
        for ci in cons
            push!(rows, ci)
            object = constraint_object(ci)
            MOI.Nonlinear.add_constraint(nlp, object.func, object.set)
        end
    end
    
    MOI.Nonlinear.set_objective(nlp, objective_function(model))
    x = all_variables(model)
    nlp_backend = MOI.Nonlinear.SparseReverseMode()
    evaluator = MOI.Nonlinear.Evaluator(nlp, nlp_backend, index.(x))
    MOI.initialize(evaluator, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(evaluator)
    #print(hessian_sparsity)

    I = [i for (i, _) in hessian_sparsity]
    J = [j for (_, j) in hessian_sparsity]
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(evaluator, V, value.(x), 1.0, dual.(rows))
    H = SparseArrays.sparse(I, J, V, length(x), length(x))
    GLMakie.spy(H)
    return Matrix(fill_off_diagonal(H))
end
H_star = compute_optimal_hessian(model)
GLMakie.spy(reverse(H_star,dims=1))
UnicodePlots.spy(H_star)