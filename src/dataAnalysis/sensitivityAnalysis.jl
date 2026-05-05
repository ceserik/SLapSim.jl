
function _collect_role!(obj, out::Vector{carParameter}, role::Symbol)
    for field in fieldnames(typeof(obj))
        value = getfield(obj, field)
        if value isa carParameter
            value.role == role && push!(out, value)
        elseif value isa AbstractVector
            for item in value
                isstructtype(typeof(item)) && !isa(item, carParameter) && _collect_role!(item, out, role)
            end
        elseif !isa(value, Function) && !isa(value, Number) && !isa(value, String) && !isa(value, Symbol) && isstructtype(typeof(value))
            _collect_role!(value, out, role)
        end
    end
end

function setParameters(car::Car, model::JuMP.Model)
    sensitivityParams = carParameter[]
    _collect_role!(car, sensitivityParams, :sensitivity)
    params = Dict{String, JuMP.VariableRef}()
    for p in sensitivityParams
        jump_param = @variable(model, set = Parameter(p.value))
        params[p.name] = jump_param
        p.value = jump_param
    end
    return params, sensitivityParams
end

function resetParameters(sensitivityParams::Vector{carParameter})
    for p in sensitivityParams
        p.value = JuMP.value(p.value)
    end
end

function setDesignVariables(car::Car, model::JuMP.Model)
    designs = carParameter[]
    _collect_role!(car, designs, :design)
    refs = Dict{String, JuMP.VariableRef}()
    for p in designs
        lb, ub = p.limits[1], p.limits[2]
        start_val = p.value isa Real ? Float64(p.value) : 0.0
        v = @variable(model, base_name = p.name, start = start_val)
        if lb < ub
            JuMP.set_lower_bound(v, lb)
            JuMP.set_upper_bound(v, ub)
        end
        refs[p.name] = v
        p.value = v
    end
    return refs, designs
end

function resolveDesignVariables(designs::Vector{carParameter})
    for p in designs
        p.value = JuMP.value(p.value)
    end
end

function sensitivityAnalysis(problem)
    model = problem.model
    params = problem.params
    results_path = problem.analysis.results_path

    DiffOpt.empty_input_sensitivities!(model)
    DiffOpt.set_reverse_objective(model, 1.0)
    try
        DiffOpt.reverse_differentiate!(model)
    catch e
        @warn "sensitivityAnalysis: DiffOpt differentiation failed (likely numerical dual sign issue), skipping." exception=e
        return nothing
    end

    obj_val = objective_value(model)
    names = collect(keys(params))
    sensitivities = [(DiffOpt.get_reverse_parameter(model, params[n]) * JuMP.value(params[n])) / obj_val for n in names]

    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel="Parameter", ylabel="Normalized sensitivity (% obj / % param)",
        title="Sensitivity Analysis",
        xticks=(1:length(names), names),
        xticklabelrotation=π/4)
    GLMakie.barplot!(ax, 1:length(names), sensitivities)
    display(GLMakie.Screen(), fig)

    mkpath(results_path)
    CairoMakie.save(joinpath(results_path, "sensitivity.svg"), fig)
    println("Sensitivity plot saved to $(joinpath(results_path, "sensitivity.svg"))")

    return Dict(zip(names, sensitivities))
end
