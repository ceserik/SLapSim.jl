# =============================================================================
# transcriptionImage.jl
#
# Produces a single static, presentation-ready image that *symbolises direct
# transcription*: the selected track with a full vehicle drawn at every
# discretisation node of the NLP mesh. Each car shows the states (pose, heading)
# and controls (steering, tyre forces) that the optimiser holds at that node,
# emphasising that every NLP iteration carries states + controls for each point
# on the track.
#
# Run with:
#   julia --project=. src/experiments/transcriptionImage.jl
# =============================================================================

using Revise
using SLapSim
using CairoMakie          # crisp, headless vector/raster export for slides

CairoMakie.activate!()

# ---------------------------------------------------------------------------
# Configuration — tweak these for the slide you want.
# ---------------------------------------------------------------------------
const N_SEGMENTS   = 20                 # NLP mesh intervals  ->  N_SEGMENTS+1 nodes (cars)
const POL_ORDER    = 3                  # collocation polynomial order per segment
const OUTPUT_STEM  = "sync/transcription"  # ".png" and ".pdf" are appended
const SHOW_RACELINE   = true            # draw the continuous optimal line under the cars
const COLOR_BY_SPEED  = true            # colour the raceline by speed (+ colourbar)
const SHOW_STATION_DOTS = true          # mark each discretisation station on the centreline
const HIDE_AXES       = true            # clean, decoration-free look for slides

# ---------------------------------------------------------------------------
# Scenario — pick a track + car. (These match a known-good configuration.)
# ---------------------------------------------------------------------------
track = doubleTurn(false, 0.1)
# Alternatives:
#   track = kml2track("tracks/FSCZ.kml", false, true)   # real Formula Student circuit
#   track = skidpad(false)
car = createTwintrack(true, track)

ipopt_attrs = Dict{String,Any}(
    "linear_solver" => "mumps",
    "alpha_for_y"   => "safer-min-dual-infeas",
    "print_level"   => 0,
)

exp = Experiment(
    car        = car,
    track      = track,
    discipline = Open(v_start = 5.0),
    solver     = IpoptBackend(performSensitivity = false, attributes = ipopt_attrs),
    mesh_refinement = MeshRefinementConfig(
        tol            = 1e-1,
        method         = :h,
        error_method   = :ode,
        max_iterations = 0,             # single solve on exactly N_SEGMENTS intervals (no refinement)
        segments       = N_SEGMENTS,
        pol_order      = POL_ORDER,
        variant        = "Radau",
    ),
    # Everything off: run headless, no pop-up windows / no animation.
    analysis = AnalysisConfig(
        plot_path = false, plot_states = false, plot_controls = false,
        plot_jacobian = false, plot_hessian = false, animate = false,
        time_simulation = false, plot_initialization = false,
    ),
)

# ---------------------------------------------------------------------------
# Geometry helper: global pose (x, y, heading) of the car at arc-length s.
# Mirrors the convention used by the animation:
#   state[3] = heading psi,  state[5] = lateral offset n from the centreline.
# ---------------------------------------------------------------------------
function _pose_at(res, track, s)
    st        = Float64.(res.states(s))
    fc        = track.fcurve(s)
    theta     = fc[2]
    cx        = fc[3] - st[5] * sin(theta)
    cy        = fc[4] + st[5] * cos(theta)
    return Float64(cx), Float64(cy), Float64(st[3]), st
end

# ---------------------------------------------------------------------------
# Render the transcription image from a solved experiment.
# ---------------------------------------------------------------------------
function render_transcription(exp; output_stem = OUTPUT_STEM,
                              show_raceline = SHOW_RACELINE,
                              color_by_speed = COLOR_BY_SPEED,
                              show_station_dots = SHOW_STATION_DOTS,
                              hide_axes = HIDE_AXES)
    track = exp.track
    car   = exp.car
    res   = exp.optiResult
    res === nothing && error("Experiment has no optiResult — solve it first with run_experiment!(exp).")

    # Discretisation stations = the NLP mesh nodes where states/controls live.
    stations = sort!(unique!(collect(Float64, res.segment_edges)))
    isempty(stations) && (stations = sort!(unique!(collect(Float64, res.path))))

    wb = car.chassis.wheelbase.value
    tw = car.chassis.track.value
    Fz_static = car.chassis.mass.value * 9.81 / length(car.wheelAssemblies)

    CairoMakie.activate!()
    fig = Figure(size = (1900, 1150), backgroundcolor = :white)
    ax  = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = :white,
               title = "Direct transcription — states & controls at every discretisation node",
               titlesize = 24)

    # (1) Track: dashed centreline + solid boundaries.
    plotTrack(track; ax = ax)

    # (2) Continuous optimal raceline (fine sampling of the interpolated solution).
    if show_raceline
        s_fine = collect(LinRange(stations[1], stations[end], 600))
        xs = similar(s_fine); ys = similar(s_fine); vs = similar(s_fine)
        for (i, s) in enumerate(s_fine)
            cx, cy, _, st = _pose_at(res, track, s)
            xs[i] = cx; ys[i] = cy; vs[i] = st[1]
        end
        if color_by_speed
            pl = lines!(ax, xs, ys; color = vs, colormap = :turbo, linewidth = 4)
            Colorbar(fig[1, 2], pl, label = "speed vₓ [m/s]", labelsize = 18)
        else
            lines!(ax, xs, ys; color = (:dodgerblue, 0.9), linewidth = 4)
        end
    end

    # (3) A full vehicle at every discretisation node (fresh plot objects per car
    #     so all of them persist simultaneously in one static image).
    cog_x = Float64[]; cog_y = Float64[]
    for s in stations
        cx, cy, psi, state = _pose_at(res, track, s)
        ctrl = Float64.(res.controls(s))

        # Populate the car so the drawing reflects this node's states/controls
        # (steering angle, tyre forces / friction ellipses, load transfer …).
        car.stateMapping(state)
        car.controlMapping(ctrl)
        try
            car.carFunction(track, nothing)
        catch err
            @warn "carFunction failed at s=$s; drawing with previous forces" exception = err
        end

        chassis_obs = car.chassis.setupObservables(ax)
        wa_obs = [wa.setupObservables(ax, car.drivetrain.tires[i])
                  for (i, wa) in enumerate(car.wheelAssemblies)]
        aero_obs = car.aero.setupObservables(ax)

        car.chassis.updateObservables(chassis_obs, cx, cy, psi)
        for (j, wa) in enumerate(car.wheelAssemblies)
            wa.updateObservables(wa_obs[j], car.drivetrain.tires[j], cx, cy, psi, Fz_static)
        end
        car.aero.updateObservables(aero_obs, cx, cy, psi, wb, tw)

        push!(cog_x, cx); push!(cog_y, cy)
    end

    # Centre-of-gravity marker on each car (the "state" sample point).
    scatter!(ax, cog_x, cog_y; color = :yellow, marker = :cross,
             markersize = 10, strokewidth = 1.0, strokecolor = :black)

    # (4) Discretisation stations on the track centreline (the s_i grid).
    if show_station_dots
        xc = [track.fcurve(s)[3] for s in stations]
        yc = [track.fcurve(s)[4] for s in stations]
        scatter!(ax, xc, yc; color = :red, marker = :circle, markersize = 8,
                 strokewidth = 0.5, strokecolor = :black)
    end

    hide_axes && (hidedecorations!(ax); hidespines!(ax))
    resize_to_layout!(fig)

    mkpath(dirname(output_stem))
    png_path = output_stem * ".png"
    pdf_path = output_stem * ".pdf"
    CairoMakie.save(png_path, fig; px_per_unit = 2)   # high-resolution raster
    CairoMakie.save(pdf_path, fig)                    # vector version for slides
    @info "Transcription image saved" nodes = length(stations) png = png_path pdf = pdf_path
    return fig
end

# ---------------------------------------------------------------------------
# Solve + render.
# ---------------------------------------------------------------------------
run_experiment!(exp)
fig = render_transcription(exp)

nothing
