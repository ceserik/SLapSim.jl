# =============================================================================
# solveAnimation.jl
#
# Produces an animation of the NLP *solving process*. The solver is paused after
# every iteration and a frame is rendered, showing a full vehicle drawn at every
# discretisation node of the mesh. As the iterations progress the fleet of cars
# morphs from the initial guess into the optimal trajectory — visualising how the
# optimiser holds states + controls for each point on the track and refines them.
#
# How the "pause every iteration" works
# --------------------------------------
# We run the optimisation ONCE and hook into Ipopt's intermediate callback
# (`Ipopt.CallbackFunction`), which fires after every solver iteration. On each
# call we read the current primal iterate, rebuild it into an interpolated
# trajectory, and stash it as an animation frame. This captures the entire solve
# in a single optimisation — no wasteful re-solving. Enabled by setting
# `exp.iterate_hook` together with `IpoptBackend(direct = true)`.
#
# Run with:
#   julia --project=. src/experiments/solveAnimation.jl
# =============================================================================

using Revise
using SLapSim
using CairoMakie
import SLapSim: solve_success

CairoMakie.activate!()

# ---------------------------------------------------------------------------
# Configuration — tweak these for the slide you want.
# ---------------------------------------------------------------------------
const N_SEGMENTS   = 15                       # mesh intervals -> N_SEGMENTS+1 nodes (cars)
const POL_ORDER    = 3                        # collocation polynomial order per segment
# Iterations change the trajectory a lot early on and barely at the end, so for
# the animation we show every iteration up to DENSE_UNTIL, then every COARSE_STRIDE
# -th one. (All iterations are captured; this only controls which become frames.)
const DENSE_UNTIL  = 15                       # show every iteration up to here …
const DENSE_STRIDE = 1
const COARSE_STRIDE = 3                       # … then every 3rd iteration afterwards
const FRAMERATE    = 10                       # animation frames per second
const HOLD_START   = 6                        # repeat first frame (initial guess) N times
const HOLD_END     = 15                       # repeat last frame (converged) N times
const OUTPUT       = "sync/animations/solve_process.mp4"
const SHOW_RACELINE = true                    # draw the current iterate's continuous line
const HIDE_AXES     = true                    # clean, decoration-free look for slides

# ---------------------------------------------------------------------------
# Scenario — same track + car as the transcription image for consistency.
# ---------------------------------------------------------------------------
track = doubleTurn(false, 0.1)
car   = createTwintrack(true, track)

ipopt_attrs = Dict{String,Any}(
    "linear_solver" => "mumps",
    "alpha_for_y"   => "safer-min-dual-infeas",
    "print_level"   => 0,
)

function make_experiment()
    Experiment(
        car        = car,
        track      = track,
        discipline = Open(v_start = 5.0),
        # direct = true is required so the per-iterate callback can read the raw
        # Ipopt iterate vector with a 1:1 variable-to-column mapping.
        solver     = IpoptBackend(performSensitivity = false, direct = true,
                                  attributes = copy(ipopt_attrs)),
        mesh_refinement = MeshRefinementConfig(
            tol            = 1e-1,
            method         = :h,
            error_method   = :ode,
            max_iterations = 0,               # single solve on exactly N_SEGMENTS intervals (no mesh refinement)
            segments       = N_SEGMENTS,
            pol_order      = POL_ORDER,
            variant        = "Radau",
        ),
        # Everything off: run headless, no pop-up windows / no animation / no plots.
        analysis = AnalysisConfig(
            plot_path = false, plot_states = false, plot_controls = false,
            plot_jacobian = false, plot_hessian = false, animate = false,
            time_simulation = false, plot_initialization = false,
        ),
    )
end

# ---------------------------------------------------------------------------
# Geometry helper: global pose (x, y, heading) of the car at arc-length s.
#   state[3] = heading psi,  state[5] = lateral offset n from the centreline,
#   state[6] = time.
# ---------------------------------------------------------------------------
function _pose_at(res, track, s)
    st    = Float64.(res.states(s))
    fc    = track.fcurve(s)
    theta = fc[2]
    cx    = fc[3] - st[5] * sin(theta)
    cy    = fc[4] + st[5] * cos(theta)
    return Float64(cx), Float64(cy), Float64(st[3]), st
end

# Return a writable path: prefer `path`, but if it exists and is locked (e.g. the
# previous video is open in a player), fall back to path_1, path_2, …
function _free_path(path)
    isfile(path) || return path
    try
        rm(path; force = true)   # deletable → reuse the canonical name
        return path
    catch
        base, ext = splitext(path)
        i = 1
        while isfile("$(base)_$(i)$(ext)")
            i += 1
        end
        newp = "$(base)_$(i)$(ext)"
        @warn "output is locked (open elsewhere?); writing to $newp instead" locked = path
        return newp
    end
end

# Select which captured iterations become frames: dense early, coarse later.
# Always keeps iteration 0 and the final (converged) iteration.
function _select_frames(all_frames)
    isempty(all_frames) && return all_frames
    lastk = all_frames[end].iter
    keep(k) = k < DENSE_UNTIL ? (k % DENSE_STRIDE == 0) : (k % COARSE_STRIDE == 0)
    sel = [f for f in all_frames if keep(f.iter) || f.iter == 0 || f.iter == lastk]
    return sel
end

# ---------------------------------------------------------------------------
# Capture phase: solve ONCE and grab every iterate through Ipopt's intermediate
# callback (wired up via exp.iterate_hook). Returns (frames, converged) where
# frames is a Vector of (iter, res, converged).
# ---------------------------------------------------------------------------
function capture_iterates()
    exp = make_experiment()
    all_frames = NamedTuple[]
    exp.iterate_hook = (it, res) -> push!(all_frames, (iter = it, res = res, converged = false))

    println("Solving once, capturing every iterate via Ipopt callback…")
    run_experiment!(exp)
    converged = solve_success(exp.model)
    println("Solve produced $(length(all_frames)) iterates" *
            (converged ? " (converged)" : " (did NOT converge)"))

    # Mark the final captured iterate as the converged one for labelling.
    if !isempty(all_frames)
        last = all_frames[end]
        all_frames[end] = (iter = last.iter, res = last.res, converged = converged)
    end

    frames = _select_frames(all_frames)
    println("Using $(length(frames)) of $(length(all_frames)) iterates as frames.")
    return frames, converged
end

# ---------------------------------------------------------------------------
# Render phase: animate the captured iterates. A car is drawn at every mesh node
# and updated in place each frame via the observable-based drawing primitives.
# ---------------------------------------------------------------------------
function render_animation(frames, track, car; output = OUTPUT, framerate = FRAMERATE,
                          show_raceline = SHOW_RACELINE, hide_axes = HIDE_AXES)
    isempty(frames) && error("No frames captured — nothing to animate.")

    wb        = car.chassis.wheelbase.value
    tw        = car.chassis.track.value
    Fz_static = car.chassis.mass.value * 9.81 / length(car.wheelAssemblies)

    # Mesh is fixed across iterations, so the node stations are the same for all frames.
    ref      = frames[end].res
    stations = sort!(unique!(collect(Float64, ref.segment_edges)))
    isempty(stations) && (stations = sort!(unique!(collect(Float64, ref.path))))

    CairoMakie.activate!()
    fig = Figure(size = (1900, 1150), backgroundcolor = :white)
    title_obs = Observable("NLP iteration 0")
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = :white,
              title = title_obs, titlesize = 26)

    # Static track (dashed centreline + solid boundaries).
    plotTrack(track; ax = ax)

    # Fixed view so nothing jumps between frames.
    s_lim = LinRange(stations[1], stations[end], 600)
    xc = [track.fcurve(s)[3] for s in s_lim]
    yc = [track.fcurve(s)[4] for s in s_lim]
    wmax = max(maximum(track.widthL), maximum(track.widthR))
    pad  = wmax + max(wb, tw) * 1.5
    xlims!(ax, minimum(xc) - pad, maximum(xc) + pad)
    ylims!(ax, minimum(yc) - pad, maximum(yc) + pad)

    # Continuous raceline for the current iterate.
    raceline = Observable(Point2f[])
    show_raceline && lines!(ax, raceline; color = (:dodgerblue, 0.9), linewidth = 3)

    # One persistent set of car observables per node (created once, updated each frame).
    chassis_obs = [car.chassis.setupObservables(ax) for _ in stations]
    wa_obs      = [[wa.setupObservables(ax, car.drivetrain.tires[i])
                    for (i, wa) in enumerate(car.wheelAssemblies)] for _ in stations]
    aero_obs    = [car.aero.setupObservables(ax) for _ in stations]

    # Centre-of-gravity markers (one point per node).
    cog_obs = Observable(Point2f[])
    scatter!(ax, cog_obs; color = :yellow, marker = :cross, markersize = 9,
             strokewidth = 1.0, strokecolor = :black)

    hide_axes && (hidedecorations!(ax); hidespines!(ax))
    resize_to_layout!(fig)

    # Update everything to show a single captured iterate.
    function draw_iterate!(frame)
        res = frame.res
        if show_raceline
            pts = Vector{Point2f}(undef, length(s_lim))
            for (i, s) in enumerate(s_lim)
                cx, cy, _, _ = _pose_at(res, track, s)
                pts[i] = Point2f(cx, cy)
            end
            raceline[] = pts
        end

        cogs = Vector{Point2f}(undef, length(stations))
        for (si, s) in enumerate(stations)
            cx, cy, psi, state = _pose_at(res, track, s)
            ctrl = Float64.(res.controls(s))
            # Populate the car so the drawing reflects this node's states/controls.
            car.stateMapping(state)
            car.controlMapping(ctrl)
            try
                car.carFunction(track, nothing)
            catch
                # Early iterates can hit non-physical intermediate points; keep previous forces.
            end
            car.chassis.updateObservables(chassis_obs[si], cx, cy, psi)
            for (j, wa) in enumerate(car.wheelAssemblies)
                wa.updateObservables(wa_obs[si][j], car.drivetrain.tires[j], cx, cy, psi, Fz_static)
            end
            car.aero.updateObservables(aero_obs[si], cx, cy, psi, wb, tw)
            cogs[si] = Point2f(cx, cy)
        end
        cog_obs[] = cogs

        T = Float64(res.states(stations[end])[6])
        stage = frame.converged ? "converged" : "NLP iteration $(frame.iter)"
        title_obs[] = "Direct transcription — $stage,  lap time T = $(round(T, digits = 3)) s"
    end

    # Frame schedule: linger on the initial guess and the converged result.
    order = vcat(fill(1, HOLD_START), collect(eachindex(frames)), fill(lastindex(frames), HOLD_END))

    mkpath(dirname(output))
    output = _free_path(output)   # avoid EBUSY if the previous video is open in a player
    n = length(order)
    CairoMakie.record(fig, output, 1:n; framerate = framerate) do fpos
        draw_iterate!(frames[order[fpos]])
        print("\rRendering frame $fpos / $n")
    end
    println("\nAnimation saved to $output")
    return fig
end

# ---------------------------------------------------------------------------
# Capture + render.
# ---------------------------------------------------------------------------
frames, converged = capture_iterates()
fig = render_animation(frames, track, car)

nothing
