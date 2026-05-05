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

"""
    compute_optimal_hessian(model) -> Matrix

Build the MOI nonlinear evaluator on `model` and return the Hessian of the
Lagrangian at the current primal / dual solution. Off-diagonals are mirrored
so the returned matrix is dense-symmetric.
"""
function compute_optimal_hessian(model::Model)
    rows = Any[]
    nlp = MOI.Nonlinear.Model()
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

    I = [i for (i, _) in hessian_sparsity]
    J = [j for (_, j) in hessian_sparsity]
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(evaluator, V, value.(x), 1.0, dual.(rows))
    H = SparseArrays.sparse(I, J, V, length(x), length(x))
    return Matrix(fill_off_diagonal(H))
end

"""
    plot_hessian_spy(model) -> H_star

Compute the Hessian of the Lagrangian on `model` and render a UnicodePlots spy
to stdout plus a GLMakie spy window. Returns the dense Hessian.
"""
function plot_hessian_spy(model::Model; savepath::Union{Nothing,String}=nothing)
    H_star = compute_optimal_hessian(model)
    println("hessian")
    println(UnicodePlots.spy(H_star))
    fig_hess = Figure()
    ax_hess = Axis(fig_hess[1, 1], title="Hessian", yreversed=true)
    spy!(ax_hess, SparseArrays.sparse(H_star))
    if savepath !== nothing
        CairoMakie.save(savepath, fig_hess)
    end
    display(GLMakie.Screen(), fig_hess)
    return H_star
end
