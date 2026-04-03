using Revise
using Interpolations

function _fmt_time(t::Float64)
    s = string(round(t; digits=3))
    d = length(s) - findlast('.', s)
    return "t = " * s * "0"^(3 - d) * " s"
end

function _unique_path(path::String)
    isfile(path) || return path
    base, ext = splitext(path)
    i = 1
    while isfile("$(base)_$i$ext")
        i += 1
    end
    return "$(base)_$i$ext"
end

# ─────────────────── Shared geometry helpers ──────────────────────────────────

"""
    _rotmat2d(θ)

2D rotation matrix for angle θ.
"""
_rotmat2d(θ::Float64) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

"""
    _rect_corners(cx, cy, w, h, θ)

Return 5×2 matrix of corners (closed polygon) for a rectangle centered at (cx,cy)
with width `w`, height `h`, rotated by `θ`.
"""
function _rect_corners(cx::Float64, cy::Float64, w::Float64, h::Float64, θ::Float64)
    hw, hh = w / 2, h / 2
    corners = [-hw -hh; hw -hh; hw hh; -hw hh; -hw -hh]
    R = _rotmat2d(θ)
    rotated = (R * corners')'
    return rotated .+ [cx cy]
end

# ─────────────────── Main car drawing function ────────────────────────────────

"""
    drawCar!(ax, x, y, ψ, car::Car, steering_angle=0.0)

Draw the full car at position (x,y) with heading ψ.
Each component draws itself via its own `draw!` method.
"""
function drawCar!(x, y, ψ, car::Car; ax=nothing)
    created = false
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1], aspect=DataAspect())
        created = true
    end

    draw!(ax, car.chassis, x, y, ψ)
    for (i, wa) in enumerate(car.wheelAssemblies)
        draw!(ax, wa, car.drivetrain.tires[i], x, y, ψ)
    end
    draw!(ax, car.aero, x, y, ψ, car.chassis.wheelbase.value, car.chassis.track.value)
    # CoG marker
    scatter!(ax, [x], [y]; color=:yellow, markersize=8, marker=:cross, strokewidth=1.5, strokecolor=:black)
    # Heading arrow
    arrows2d!(ax, [x], [y], [0.4 * cos(ψ)], [0.4 * sin(ψ)]; color=:orange, shaftwidth=2, tipwidth=10, tiplength=10)

    if created
        display(GLMakie.Screen(), fig)
        return fig, ax
    end
    return ax
end

# ─────────────────── Animation function ───────────────────────────────────────

"""
    _rect_points(cx, cy, w, h, θ)

Same as `_rect_corners` but returns `Vector{Point2f}` for Observable updates.
"""
function _rect_points(cx::Float64, cy::Float64, w::Float64, h::Float64, θ::Float64)
    c = _rect_corners(cx, cy, w, h, θ)
    return Point2f.(eachrow(c))
end

"""
    _setup_panel(ax, track, car, follow_car, view_radius, xc, yc, thc, xl, yl, xr, yr)

Set up one animation panel (axis with car observables). Returns a NamedTuple of observables.
"""
function _setup_panel(ax::Axis, track::Track, car::Car, follow_car::Bool, view_radius::Float64, xc::Vector{Float64}, yc::Vector{Float64}, thc::Vector{Float64}, xl::Vector{Float64}, yl::Vector{Float64}, xr::Vector{Float64}, yr::Vector{Float64}; cam_offset::Float64=2.0)
    track_center_obs = nothing
    track_left_obs = nothing
    track_right_obs = nothing

    # Pre-allocate buffers for follow-cam track transforms
    n_track = length(xc)
    buf_center = Vector{Point2f}(undef, n_track)
    buf_left   = Vector{Point2f}(undef, n_track)
    buf_right  = Vector{Point2f}(undef, n_track)

    if follow_car
        track_center_obs = Observable(Point2f.(xc, yc))
        track_left_obs   = Observable(Point2f.(xl, yl))
        track_right_obs  = Observable(Point2f.(xr, yr))
        lines!(ax, track_center_obs; linestyle=:dash, linewidth=1)
        lines!(ax, track_left_obs; color=:black, linewidth=1)
        lines!(ax, track_right_obs; color=:black, linewidth=1)
        xlims!(ax, -view_radius, view_radius)
        ylims!(ax, -view_radius + cam_offset, view_radius + cam_offset)
    else
        plotTrack(track, b_plotStartEnd=false, ax=ax)
    end

    # Each component sets up its own observables
    chassis_obs = car.chassis.setupObservables(ax)
    wa_obs = [wa.setupObservables(ax, car.drivetrain.tires[i]) for (i, wa) in enumerate(car.wheelAssemblies)]
    aero_obs = car.aero.setupObservables(ax)

    cog_obs = Observable([Point2f(0, 0)])
    scatter!(ax, cog_obs; color=:yellow, markersize=8, marker=:cross, strokewidth=1.5, strokecolor=:black)

    trail_obs = Observable(Point2f[])
    lines!(ax, trail_obs; color=(:dodgerblue, 0.5), linewidth=2)

    return (chassis=chassis_obs, wheel_assemblies=wa_obs, aero=aero_obs,
            cog=cog_obs, trail=trail_obs,
            track_center=track_center_obs, track_left=track_left_obs, track_right=track_right_obs,
            trail_global=Point2f[], follow=follow_car,
            buf_center=buf_center, buf_left=buf_left, buf_right=buf_right)
end

"""
    _update_panel!(panel, car, cx, cy, ψ, xc, yc, xl, yl, xr, yr, view_radius)

Update one panel's observables for the current frame.
Reads steering angles and forces directly from the car struct (set by controlMapping/carFunction).
"""
function _update_panel!(panel::NamedTuple, car::Car, cx::Float64, cy::Float64, ψ::Float64, xc::Vector{Float64}, yc::Vector{Float64}, xl::Vector{Float64}, yl::Vector{Float64}, xr::Vector{Float64}, yr::Vector{Float64}, view_radius::Float64)
    wb = car.chassis.wheelbase.value
    tw = car.chassis.track.value
    Fz_static = car.chassis.mass.value * 9.81 / length(car.wheelAssemblies)

    if panel.follow
        Rcam = _rotmat2d(π/2 - ψ)
        r11, r12, r21, r22 = Rcam[1,1], Rcam[1,2], Rcam[2,1], Rcam[2,2]

        @inbounds for k in eachindex(xc)
            dx, dy = xc[k] - cx, yc[k] - cy
            panel.buf_center[k] = Point2f(r11*dx + r12*dy, r21*dx + r22*dy)
        end
        @inbounds for k in eachindex(xl)
            dx, dy = xl[k] - cx, yl[k] - cy
            panel.buf_left[k] = Point2f(r11*dx + r12*dy, r21*dx + r22*dy)
        end
        @inbounds for k in eachindex(xr)
            dx, dy = xr[k] - cx, yr[k] - cy
            panel.buf_right[k] = Point2f(r11*dx + r12*dy, r21*dx + r22*dy)
        end
        panel.track_center[] = panel.buf_center
        panel.track_left[]   = panel.buf_left
        panel.track_right[]  = panel.buf_right

        car_ψ = π/2
        car.chassis.updateObservables(panel.chassis, 0.0, 0.0, car_ψ)
        for (j, wa) in enumerate(car.wheelAssemblies)
            wa.updateObservables(panel.wheel_assemblies[j], car.drivetrain.tires[j], 0.0, 0.0, car_ψ, Fz_static)
        end
        car.aero.updateObservables(panel.aero, 0.0, 0.0, car_ψ, wb, tw)

        panel.cog[] = [Point2f(0, 0)]

        push!(panel.trail_global, Point2f(cx, cy))
        trail_buf = panel.trail[]
        resize!(trail_buf, length(panel.trail_global))
        @inbounds for k in eachindex(panel.trail_global)
            dx, dy = panel.trail_global[k][1] - cx, panel.trail_global[k][2] - cy
            trail_buf[k] = Point2f(r11*dx + r12*dy, r21*dx + r22*dy)
        end
        notify(panel.trail)
    else
        car.chassis.updateObservables(panel.chassis, cx, cy, ψ)
        for (j, wa) in enumerate(car.wheelAssemblies)
            wa.updateObservables(panel.wheel_assemblies[j], car.drivetrain.tires[j], cx, cy, ψ, Fz_static)
        end
        car.aero.updateObservables(panel.aero, cx, cy, ψ, wb, tw)

        panel.cog[] = [Point2f(cx, cy)]

        push!(panel.trail[], Point2f(cx, cy))
        notify(panel.trail)
    end
end

# ── Precompute track boundaries (shared by all animation functions) ──
function _precompute_track(track::Track)
    s_track = track.sampleDistances
    vals = track.fcurve.(s_track)
    xc  = getindex.(vals, 3)
    yc  = getindex.(vals, 4)
    thc = getindex.(vals, 2)
    xl = xc .- track.widthL .* sin.(thc)
    yl = yc .+ track.widthL .* cos.(thc)
    xr = xc .+ track.widthR .* sin.(thc)
    yr = yc .- track.widthR .* cos.(thc)
    return xc, yc, thc, xl, yl, xr, yr
end

# ── Animation loop (shared) ──
function _animate_loop!(fig::Figure, panels::Vector{<:NamedTuple}, track::Track, result, car::Car, speedup::Float64, time_obs::Observable{String}, xc::Vector{Float64}, yc::Vector{Float64}, xl::Vector{Float64}, yl::Vector{Float64}, xr::Vector{Float64}, yr::Vector{Float64}, view_radius::Float64)
    is_matrix = result.states isa AbstractMatrix
    n_frames = is_matrix ? size(result.states, 1) : length(result.path)

    for i in 1:n_frames
        if is_matrix
            state = result.states[i, :]
            ctrl = result.controls[min(i, size(result.controls, 1)), :]
            s = result.path[i]
        else
            s = result.path[i]
            state = result.states(s)
            ctrl = result.controls(s)
        end

        # Update car struct — sets all internal parameters, forces, steering, etc.
        car.stateMapping(state)
        car.controlMapping(ctrl)
        car.carFunction(track, nothing)

        ψ = state[3]
        n = state[5]

        fc = track.fcurve(s)
        theta = fc[2]
        cx = fc[3] - n * sin(theta)
        cy = fc[4] + n * cos(theta)

        for panel in panels
            _update_panel!(panel, car, cx, cy, ψ, xc, yc, xl, yl, xr, yr, view_radius)
        end

        time_obs[] = _fmt_time(state[6])
    end
end

"""
    animateCar(track::Track, result, car::Car; fps=30, speedup=1.0, follow_car=false, view_radius=10.0)

Animate the car driving along the track.
- `follow_car=false` (default): global frame, full track visible.
- `follow_car=true`: camera follows the car, showing a window of `view_radius` around it.
"""
function animateCar(track::Track, result, car::Car; fps=30, speedup=1.0, follow_car=false, view_radius=10.0, cam_offset=0.0)
    xc, yc, thc, xl, yl, xr, yr = _precompute_track(track)

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())
    panel = _setup_panel(ax, track, car, follow_car, view_radius, xc, yc, thc, xl, yl, xr, yr; cam_offset)

    time_obs = Observable("t = 0.000 s")
    timer_pos = Point2f(follow_car ? view_radius * 0.95 : maximum(track.x),
                        follow_car ? view_radius * 0.9  : maximum(track.y))
    text!(ax, timer_pos; text=time_obs, fontsize=16, align=(:right, :top))

    colsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)
    display(GLMakie.Screen(), fig)

    _animate_loop!(fig, [panel], track, result, car, speedup, time_obs, xc, yc, xl, yl, xr, yr, view_radius)
    return fig
end

"""
    animateCarDual(track::Track, result, car::Car; fps=30, speedup=1.0, view_radius=10.0)

Animate with two views side by side: global frame (left) and follow-car frame (right).
"""
function animateCarDual(track::Track, result, car::Car; fps=30, speedup=1.0, view_radius=10.0, savepath=nothing, cam_offset=0.0)
    println("starting animation")
    xc, yc, thc, xl, yl, xr, yr = _precompute_track(track)
    #fig = Figure(size=(3840, 2160))
    fig = Figure(size=(1920, 1080))
    ax_global = Axis(fig[1, 1], aspect=DataAspect(), title="Global")
    ax_follow = Axis(fig[1, 2], aspect=DataAspect(), title="Car Frame")

    panel_global = _setup_panel(ax_global, track, car, false, view_radius, xc, yc, thc, xl, yl, xr, yr)
    panel_follow = _setup_panel(ax_follow, track, car, true,  view_radius, xc, yc, thc, xl, yl, xr, yr; cam_offset)

    time_obs = Observable("t = 0.000 s")
    Label(fig[2, 1:2], time_obs; fontsize=20, halign=:center)
    rowsize!(fig.layout, 2, Fixed(30))

    panels = [panel_global, panel_follow]

    if savepath !== nothing
        savepath = _unique_path(savepath)

        # Build uniform time grid for smooth playback
        s_path = result.path
        times = [result.states isa AbstractMatrix ? result.states[i, 6] : result.states(s_path[i])[6] for i in eachindex(s_path)]
        times = accumulate(max, times)  # enforce monotonicity
        unique_mask = [true; diff(times) .> 0]
        t_unique = times[unique_mask]
        s_unique = s_path[unique_mask]
        time2s = extrapolate(interpolate((t_unique,), s_unique, Gridded(Linear())), Flat())

        t_start, t_end = t_unique[1], t_unique[end]
        dt = 1.0 / fps
        t_uniform = collect(t_start:dt/speedup:t_end)

        n_total = length(t_uniform)
        record(fig, savepath; framerate=fps, compression=5) do io
            for (i, t) in enumerate(t_uniform)
                s = time2s(t)
                state = result.states(s)
                ctrl = result.controls(s)

                car.stateMapping(state)
                car.controlMapping(ctrl)
                car.carFunction(track, nothing)

                ψ = state[3]
                n_val = state[5]
                fc = track.fcurve(s)
                theta = fc[2]
                cx = fc[3] - n_val * sin(theta)
                cy = fc[4] + n_val * cos(theta)

                for panel in panels
                    _update_panel!(panel, car, cx, cy, ψ, xc, yc, xl, yl, xr, yr, view_radius)
                end
                time_obs[] = _fmt_time(t)
                recordframe!(io)
                print("\rRendering: $i / $n_total frames")
            end
        end
        println()
        println("Animation saved to $savepath")
    else
        display(GLMakie.Screen(), fig)
        _animate_loop!(fig, panels, track, result, car, speedup, time_obs, xc, yc, xl, yl, xr, yr, view_radius)
    end
    return fig
end

function _update_frame!(panels::Vector{<:NamedTuple}, track::Track, result, car::Car, time_obs::Observable{String}, i::Int, is_matrix::Bool, xc::Vector{Float64}, yc::Vector{Float64}, xl::Vector{Float64}, yl::Vector{Float64}, xr::Vector{Float64}, yr::Vector{Float64}, view_radius::Float64)
    if is_matrix
        state = result.states[i, :]
        ctrl = result.controls[min(i, size(result.controls, 1)), :]
        s = result.path[i]
    else
        s = result.path[i]
        state = result.states(s)
        ctrl = result.controls(s)
    end

    car.stateMapping(state)
    car.controlMapping(ctrl)
    car.carFunction(track, nothing)

    ψ = state[3]
    n = state[5]
    fc = track.fcurve(s)
    theta = fc[2]
    cx = fc[3] - n * sin(theta)
    cy = fc[4] + n * cos(theta)

    for panel in panels
        _update_panel!(panel, car, cx, cy, ψ, xc, yc, xl, yl, xr, yr, view_radius)
    end
    time_obs[] = _fmt_time(state[6])
end
