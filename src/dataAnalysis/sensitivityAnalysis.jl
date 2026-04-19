
function _collect_tunable!(obj, tunables::Vector{carParameter})
    for field in fieldnames(typeof(obj))
        value = getfield(obj, field)
        if value isa carParameter
            value.role == :tunable && push!(tunables, value)
        elseif value isa AbstractVector
            for item in value
                isstructtype(typeof(item)) && !isa(item, carParameter) && _collect_tunable!(item, tunables)
            end
        elseif !isa(value, Function) && !isa(value, Number) && !isa(value, String) && !isa(value, Symbol) && isstructtype(typeof(value))
            _collect_tunable!(value, tunables)
        end
    end
end

function setParameters(car::Car, model::JuMP.Model)
    tunables = carParameter[]
    _collect_tunable!(car, tunables)
    params = Dict{String, JuMP.VariableRef}()
    for p in tunables
        jump_param = @variable(model, set = Parameter(p.value))
        params[p.name] = jump_param
        p.value = jump_param
    end
    return params, tunables
end

function resetParameters(tunables::Vector{carParameter})
    for p in tunables
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
