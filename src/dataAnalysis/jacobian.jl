using JuMP
import Ipopt
import LinearAlgebra
import Random
import SparseArrays
import MathOptInterface as MOI


function compute_optimal_jacobian(model::Model)
    rows = Any[]
    nlp = MOI.Nonlinear.Model()
    for (F, S) in list_of_constraint_types(model)
        for ci in all_constraints(model, F, S)
            if !(F <: VariableRef)
                push!(rows, ci)
                object = constraint_object(ci)
                MOI.Nonlinear.add_constraint(nlp, object.func, object.set)
            end
        end
    end
    MOI.Nonlinear.set_objective(nlp, objective_function(model))
    x = all_variables(model)
    backend = MOI.Nonlinear.SparseReverseMode()
    evaluator = MOI.Nonlinear.Evaluator(nlp, backend, index.(x))
    # Initialize the Jacobian
    MOI.initialize(evaluator, [:Jac])
    # Query the Jacobian structure
    sparsity = MOI.jacobian_structure(evaluator)
    I, J, V = first.(sparsity), last.(sparsity), zeros(length(sparsity))
    # Query the Jacobian values
    MOI.eval_constraint_jacobian(evaluator, V, value.(x))
    return SparseArrays.sparse(I, J, V, length(rows), length(x))
end

jacobian = compute_optimal_jacobian(model)
GLMakie.spy(reverse(jacobian,dims=1))
UnicodePlots.spy(jacobian)