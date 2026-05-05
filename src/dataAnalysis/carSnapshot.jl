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
        elseif val isa AbstractVector{VarEntry}
            for entry in val
                p = first(entry.targets).first
                p isa carParameter && (result[entry.name] = p)
            end
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
function snapshot_car(car::Car, optiResult, track::Track; sample_step=0.3)
    sample_points = isnothing(sample_step) ? optiResult.path : collect(optiResult.path[1]:sample_step:optiResult.path[end])
    snapshots = CarSnapshot[]
    for s in sample_points
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

# --- Role-based VarEntry accessors ---

"""
    get_entries_by_role(car, role) -> Vector{VarEntry}

Return all VarEntries across state and control descriptors whose role matches.
"""
function get_entries_by_role(car::Car, role::Symbol)
    all_entries = vcat(car.carParameters.state_descriptor, car.carParameters.control_descriptor)
    return filter(e -> e.role == role, all_entries)
end

get_states(car::Car)   = get_entries_by_role(car, :state)
get_controls(car::Car) = get_entries_by_role(car, :control)

"""
    _varentry_to_key(entry, params) -> String or Pair{String,Int}

Map a VarEntry to the snapshot parameter key used by plot_parameters.
"""
function _varentry_to_key(entry::VarEntry, params::Dict{String, carParameter})
    target_param, idx = first(entry.targets)
    for (k, v) in params
        if v === target_param
            return idx == 0 ? k : k => idx
        end
    end
    return entry.name
end

"""
    plot_states_controls(car, optiResult, track; sample_step=0.3)

Plot all states and controls in two separate windows, using plot_parameters.
"""
function plot_states_controls(car::Car, optiResult, track::Track; sample_step=0.3)
    snapshots = snapshot_car(car, optiResult, track; sample_step=sample_step)
    snapshots_nodes = snapshot_car(car, optiResult, track; sample_step=nothing)
    params = collect_carParameters(car)
    state_keys   = [_varentry_to_key(e, params) for e in get_states(car)]
    control_keys = [_varentry_to_key(e, params) for e in get_controls(car)]

    fig_states   = plot_parameters(snapshots, car, state_keys...)
    fig_controls = plot_parameters(snapshots, car, control_keys...)

    # Overlay discrete collocation nodes as X markers
    edges = Set(optiResult.segment_edges)
    _overlay_nodes!(fig_states, snapshots_nodes, state_keys, edges)
    _overlay_nodes!(fig_controls, snapshots_nodes, control_keys, edges)

    return fig_states, fig_controls
end

function _overlay_nodes!(fig, snapshots_nodes, keys, edges::Set)
    for (i, key) in enumerate(keys)
        ax = contents(fig[i, 1])[1]
        param_key = key isa Pair ? key[1] : key
        idx = key isa Pair ? key[2] : 0
        s, vals = get_scalar_values(snapshots_nodes, param_key, idx)
        is_edge = [si ∈ edges for si in s]
        # Segment edges in red, interior collocation points in black
        if any(is_edge)
            scatter!(ax, s[is_edge], vals[is_edge], marker=:xcross, markersize=20, color=:red)
        end
        if any(.!is_edge)
            scatter!(ax, s[.!is_edge], vals[.!is_edge], marker=:xcross, markersize=20, color=:black)
        end
    end
end

"""
    plot_parameters(snapshots, car, keys...; fig=nothing)

Plot parameters vs path distance. Each argument creates one subplot row.
To overlay multiple parameters on the same axis, pass them as a Vector.

Examples:
    # Each on its own axis:
    plot_parameters(snapshots, car, "drivetrain.tires[1].slipAngle", "drivetrain.tires[2].slipAngle")

    # Grouped on same axis:
    plot_parameters(snapshots, car, ["drivetrain.tires[1].forces" => 1, "drivetrain.tires[2].forces" => 1])

    # Mixed:
    plot_parameters(snapshots, car,
        ["drivetrain.tires[3].forces" => 1, "drivetrain.tires[4].forces" => 1],
        "drivetrain.tires[1].slipAngle",
    )
"""
function plot_parameters(snapshots::Vector{CarSnapshot}, car::Car, keys...; fig=nothing)
    if fig === nothing
        fig = Figure()
    end

    params = collect_carParameters(car)

    # Normalize: wrap bare keys into single-element vectors
    groups = [k isa AbstractVector ? k : [k] for k in keys]

    nplots = length(groups)
    axes = [Axis(fig[i, 1], ylabel=_make_ylabel(first(groups[i]), params)) for i in 1:nplots]
    linkxaxes!(axes...)
    axes[end].xlabel = "path [m]"

    for (i, group) in enumerate(groups)
        for key in group
            _plot_key!(axes[i], snapshots, params, key)
        end
        axislegend(axes[i], position=:rt)
    end

    display(GLMakie.Screen(), fig)
    return fig
end

function _plot_key!(ax, snapshots, params, key)
    param_key = key isa Pair ? key[1] : key
    prefix = _component_prefix(param_key)

    if key isa Pair
        _, idx = key
        s, vals = get_scalar_values(snapshots, param_key, idx)
        p = get(params, param_key, nothing)
        label = p !== nothing ? "$(prefix)$(p.name) $(_xyz(idx)) [$(p.unit)]" : "$(param_key) $(_xyz(idx))"
        lines!(ax, s, vals, label=label, linewidth=4)
    else
        s, raw = get_values(snapshots, key)
        if first(raw) isa AbstractVector
            p = get(params, key, nothing)
            pname = p !== nothing ? p.name : key
            punit = p !== nothing ? " [$(p.unit)]" : ""
            for j in eachindex(first(raw))
                vals_j = Float64[r[j] for r in raw]
                lines!(ax, s, vals_j, label="$(prefix)$(pname) $(_xyz(j))$(punit)", linewidth=4)
            end
        else
            vals = Float64.(raw)
            p = get(params, key, nothing)
            label = p !== nothing ? "$(prefix)$(p.name) [$(p.unit)]" : key
            lines!(ax, s, vals, label=label, linewidth=4)
        end
    end
end

"""Extract the parent component identifier from a key path, e.g.
"drivetrain.tires[3].forces" -> "tires[3] "
"drivetrain.motors[1].torque" -> "motors[1] "
"""
function _component_prefix(key::String)
    parts = split(key, ".")
    # Find the last part that contains an index like [3]
    for i in length(parts)-1:-1:1
        if occursin(r"\[\d+\]", parts[i])
            return string(parts[i]) * " "
        end
    end
    return ""
end

const _XYZ_LABELS = ["X", "Y", "Z"]
_xyz(i::Int) = i <= length(_XYZ_LABELS) ? _XYZ_LABELS[i] : "[$i]"

function _make_ylabel(key, params)
    param_key = key isa Pair ? key[1] : key
    p = get(params, param_key, nothing)
    if p !== nothing
        return "$(p.name) [$(p.unit)]"
    else
        return string(param_key)
    end
end
