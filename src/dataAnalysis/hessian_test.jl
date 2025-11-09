using JuMP
import Ipopt
import LinearAlgebra
import Random
import SparseArrays
import MathOptInterface as MOI


#model = Model(Ipopt.Optimizer)
#set_silent(model)
#@variable(model, x[i=1:2], start = -i)
#@constraint(model, g_1, x[1]^2 <= 1)
#@constraint(model, g_2, (x[1] + x[2])^2 <= 2)
#@objective(model, Min, (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2)
#optimize!(model)
#assert_is_solved_and_feasible(model)



rows = Any[]
nlp = MOI.Nonlinear.Model()
for (F, S) in list_of_constraint_types(model)
    if F <: VariableRef
        continue  # Skip variable bounds
    end
    for ci in all_constraints(model, F, S)
        push!(rows, ci)
        object = constraint_object(ci)
        MOI.Nonlinear.add_constraint(nlp, object.func, object.set)
    end
end
MOI.Nonlinear.set_objective(nlp, objective_function(model))
nlp

evaluator = MOI.Nonlinear.Evaluator(
    nlp,
    MOI.Nonlinear.SparseReverseMode(),
    index.(all_variables(model)),
)

MOI.initialize(evaluator, [:Hess])
hessian_sparsity = MOI.hessian_lagrangian_structure(evaluator)

I = [i for (i, _) in hessian_sparsity]
J = [j for (_, j) in hessian_sparsity]
V = zeros(length(hessian_sparsity))
n = num_variables(model)
H = SparseArrays.sparse(I, J, V, n, n)


MOI.eval_hessian_lagrangian(evaluator, V, ones(n), 1.0, ones(length(rows)))
H = SparseArrays.sparse(I, J, V, n, n)