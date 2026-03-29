using Revise
using SLapSim

"""
    collect_carParameters(obj; prefix="")

Recursively walk any object and collect all `carParameter` fields.
Returns a Dict{String, carParameter} where keys are dot-separated paths like
"drivetrain.tires[1].forces".
"""
function collect_carParameters(obj; prefix="")
    result = Dict{String, carParameter}()
    _walk!(result, obj, prefix)
    return result
end

function _walk!(result, obj, prefix)
    T = typeof(obj)
    if obj isa carParameter
        result[prefix] = obj
        return
    end
    if !isstructtype(T) || T <: AbstractString || T <: Number
        return
    end
    for fname in fieldnames(T)
        val = try getfield(obj, fname) catch; continue end
        key = isempty(prefix) ? string(fname) : prefix * "." * string(fname)
        if val isa carParameter
            result[key] = val
        elseif val isa AbstractVector
            for (i, el) in enumerate(val)
                _walk!(result, el, key * "[$i]")
            end
        elseif isstructtype(typeof(val)) && !(val isa Function) && !(val isa AbstractString) && !(val isa Number)
            _walk!(result, val, key)
        end
    end
end

"""
    CarSnapshot

Stores a deep copy of all `carParameter` values at one discretization point.
"""
struct CarSnapshot
    path::Float64
    values::Dict{String, Any}  # parameter path => value (deep copied)
end

"""
    snapshot_car(car, optiResult) -> Vector{CarSnapshot}

For each discretization point in `optiResult.path`, apply the state and control
mappings to the car, run the car function to populate all internal parameters,
and store a snapshot of every `carParameter` value.

Usage:
    snapshots = snapshot_car(car, optiResult_interp, track)
"""
function snapshot_car(car::Car, optiResult, track::Track)
    snapshots = CarSnapshot[]
    for s in optiResult.path
        x = optiResult.states(s)
        u = optiResult.controls(s)
        car.controlMapping(u)
        car.stateMapping(x)
        car.carFunction(track, nothing)

        params = collect_carParameters(car)
        vals = Dict{String, Any}()
        for (k, p) in params
            vals[k] = deepcopy(p.value)
        end
        push!(snapshots, CarSnapshot(s, vals))
    end
    return snapshots
end

"""
    get_values(snapshots, key) -> (path::Vector{Float64}, values::Vector)

Extract a single parameter across all snapshots.
`key` is the dot-separated path, e.g. "drivetrain.tires[1].forces".

For vector-valued parameters, returns a Vector{Vector}.
For scalar parameters, returns a Vector{Float64}.
"""
function get_values(snapshots::Vector{CarSnapshot}, key::String)
    path = [s.path for s in snapshots]
    vals = [s.values[key] for s in snapshots]
    return path, vals
end

"""
    get_scalar_values(snapshots, key, [index]) -> (path, values::Vector{Float64})

Extract a scalar from all snapshots. For vector parameters, provide the index.
"""
function get_scalar_values(snapshots::Vector{CarSnapshot}, key::String, index::Int=0)
    path = [s.path for s in snapshots]
    if index > 0
        vals = Float64[s.values[key][index] for s in snapshots]
    else
        vals = Float64[s.values[key] for s in snapshots]
    end
    return path, vals
end

"""
    list_parameters(snapshots) -> Vector{String}

List all available parameter keys from the first snapshot.
"""
function list_parameters(snapshots::Vector{CarSnapshot})
    return sort(collect(keys(snapshots[1].values)))
end

"""
    plot_parameters(snapshots, car, keys...; fig=nothing)

Plot one or more parameters vs path distance. Uses the car's carParameter
`.name` and `.unit` fields for axis labels and legends.

Each key is a string like "drivetrain.tires[1].forces" or a Pair like
"drivetrain.tires[1].forces" => 1 to select a vector index.

Examples:
    plot_parameters(snapshots, car, "drivetrain.tires[1].slipAngle")
    plot_parameters(snapshots, car, "drivetrain.tires[1].forces" => 1, "drivetrain.tires[2].forces" => 1)
"""
function plot_parameters(snapshots::Vector{CarSnapshot}, car::Car, keys...; fig=nothing)
    if fig === nothing
        fig = Figure()
    end

    params = collect_carParameters(car)

    nplots = length(keys)
    axes = [Axis(fig[i, 1], ylabel=_make_ylabel(keys[i], params)) for i in 1:nplots]
    linkxaxes!(axes...)
    axes[end].xlabel = "path [m]"

    for (i, key) in enumerate(keys)
        if key isa Pair
            param_key, idx = key
            s, vals = get_scalar_values(snapshots, param_key, idx)
            p = get(params, param_key, nothing)
            label = p !== nothing ? "$(p.name)[$idx]" : "$(param_key)[$idx]"
        else
            s, raw = get_values(snapshots, key)
            if first(raw) isa AbstractVector
                for j in eachindex(first(raw))
                    vals_j = Float64[r[j] for r in raw]
                    lines!(axes[i], s, vals_j, label="[$j]")
                end
                axislegend(axes[i], position=:rt)
                continue
            else
                vals = Float64.(raw)
                p = get(params, key, nothing)
                label = p !== nothing ? p.name : key
            end
        end
        lines!(axes[i], s, vals, label=label)
        axislegend(axes[i], position=:rt)
    end

    display(GLMakie.Screen(), fig)
    return fig
end

function _make_ylabel(key, params)
    param_key = key isa Pair ? key[1] : key
    p = get(params, param_key, nothing)
    if p !== nothing
        return "$(p.name) [$(p.unit)]"
    else
        return string(param_key)
    end
end
