using JuMP
import Ipopt
import LinearAlgebra
import Random
import SparseArrays
import MathOptInterface as MOI
using UnicodePlots

"""
    compute_optimal_jacobian(model) -> SparseMatrixCSC

Build the MOI nonlinear evaluator on `model` and return the constraint Jacobian
evaluated at the current primal solution.
"""
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
    MOI.initialize(evaluator, [:Jac])
    sparsity = MOI.jacobian_structure(evaluator)
    I, J, V = first.(sparsity), last.(sparsity), zeros(length(sparsity))
    MOI.eval_constraint_jacobian(evaluator, V, value.(x))
    return SparseArrays.sparse(I, J, V, length(rows), length(x))
end

"""
    plot_jacobian_spy(model) -> jacobian

Compute the Jacobian on `model` and render a UnicodePlots spy to stdout plus a
GLMakie spy window. Returns the sparse Jacobian.
"""
function plot_jacobian_spy(model::Model; savepath::Union{Nothing,String}=nothing)
    jacobian = compute_optimal_jacobian(model)
    println("jacobian")
    println(UnicodePlots.spy(jacobian))
    fig_jac = Figure()
    ax_jac = Axis(fig_jac[1, 1], title="Jacobian")
    spy!(ax_jac, SparseArrays.sparse(rotr90(jacobian)))
    if savepath !== nothing
        CairoMakie.save(savepath, fig_jac)
    end
    display(GLMakie.Screen(), fig_jac)
    return jacobian
end
