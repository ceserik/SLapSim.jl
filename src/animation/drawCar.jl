using Revise

# ─────────────────── Shared geometry helpers ──────────────────────────────────

"""
    _rotmat2d(θ)

2D rotation matrix for angle θ.
"""
_rotmat2d(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

"""
    _rect_corners(cx, cy, w, h, θ)

Return 5×2 matrix of corners (closed polygon) for a rectangle centered at (cx,cy)
with width `w`, height `h`, rotated by `θ`.
"""
function _rect_corners(cx, cy, w, h, θ)
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
function _rect_points(cx, cy, w, h, θ)
    c = _rect_corners(cx, cy, w, h, θ)
    return Point2f.(eachrow(c))
end

"""
    _setup_panel(ax, track, car, follow_car, view_radius, xc, yc, thc, xl, yl, xr, yr)

Set up one animation panel (axis with car observables). Returns a NamedTuple of observables.
"""
function _setup_panel(ax, track, car, follow_car, view_radius, xc, yc, thc, xl, yl, xr, yr)
    track_center_obs = nothing
    track_left_obs = nothing
    track_right_obs = nothing

    if follow_car
        track_center_obs = Observable(Point2f.(xc, yc))
        track_left_obs   = Observable(Point2f.(xl, yl))
        track_right_obs  = Observable(Point2f.(xr, yr))
        lines!(ax, track_center_obs; linestyle=:dash, linewidth=1)
        lines!(ax, track_left_obs; color=:black, linewidth=1)
        lines!(ax, track_right_obs; color=:black, linewidth=1)
        xlims!(ax, -view_radius, view_radius)
        ylims!(ax, -view_radius, view_radius)
    else
        plotTrack(track, b_plotStartEnd=false, ax=ax)
    end

    # Each component sets up its own observables
    chassis_obs = setup_observables!(ax, car.chassis)
    wa_obs = [setup_observables!(ax, wa, car.drivetrain.tires[i]) for (i, wa) in enumerate(car.wheelAssemblies)]
    aero_obs = setup_observables!(ax, car.aero)

    cog_obs = Observable([Point2f(0, 0)])
    scatter!(ax, cog_obs; color=:yellow, markersize=8, marker=:cross, strokewidth=1.5, strokecolor=:black)

    trail_obs = Observable(Point2f[])
    lines!(ax, trail_obs; color=(:dodgerblue, 0.5), linewidth=2)

    return (chassis=chassis_obs, wheel_assemblies=wa_obs, aero=aero_obs,
            cog=cog_obs, trail=trail_obs,
            track_center=track_center_obs, track_left=track_left_obs, track_right=track_right_obs,
            trail_global=Point2f[], follow=follow_car)
end

"""
    _update_panel!(panel, car, cx, cy, ψ, xc, yc, xl, yl, xr, yr, view_radius)

Update one panel's observables for the current frame.
Reads steering angles and forces directly from the car struct (set by controlMapping/carFunction).
"""
function _update_panel!(panel, car, cx, cy, ψ, xc, yc, xl, yl, xr, yr, view_radius)
    wb = car.chassis.wheelbase.value
    tw = car.chassis.track.value

    if panel.follow
        Rcam = _rotmat2d(π/2 - ψ)
        car_pos = [cx, cy]

        panel.track_center[] = [Point2f(Rcam * ([xc[k], yc[k]] - car_pos)) for k in eachindex(xc)]
        panel.track_left[]   = [Point2f(Rcam * ([xl[k], yl[k]] - car_pos)) for k in eachindex(xl)]
        panel.track_right[]  = [Point2f(Rcam * ([xr[k], yr[k]] - car_pos)) for k in eachindex(xr)]

        car_ψ = π/2
        update_observables!(panel.chassis, car.chassis, 0.0, 0.0, car_ψ)
        for (j, wa) in enumerate(car.wheelAssemblies)
            update_observables!(panel.wheel_assemblies[j], wa, car.drivetrain.tires[j], 0.0, 0.0, car_ψ)
        end
        update_observables!(panel.aero, car.aero, 0.0, 0.0, car_ψ, wb, tw)

        panel.cog[] = [Point2f(0, 0)]

        push!(panel.trail_global, Point2f(cx, cy))
        panel.trail[] = [Point2f(Rcam * ([p[1], p[2]] - car_pos)) for p in panel.trail_global]
    else
        update_observables!(panel.chassis, car.chassis, cx, cy, ψ)
        for (j, wa) in enumerate(car.wheelAssemblies)
            update_observables!(panel.wheel_assemblies[j], wa, car.drivetrain.tires[j], cx, cy, ψ)
        end
        update_observables!(panel.aero, car.aero, cx, cy, ψ, wb, tw)

        panel.cog[] = [Point2f(cx, cy)]

        push!(panel.trail[], Point2f(cx, cy))
        notify(panel.trail)
    end
end

# ── Precompute track boundaries (shared by all animation functions) ──
function _precompute_track(track)
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
function _animate_loop!(fig, panels, track, result, car, speedup, time_obs, xc, yc, xl, yl, xr, yr, view_radius)
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

        time_obs[] = "t = $(round(state[6]; digits=3)) s"

        if i < n_frames
            t_current = state[6]
            t_next = is_matrix ? result.states[i + 1, 6] : result.states(result.path[i + 1])[6]
            sleep(max(0.0, (t_next - t_current) / speedup))
        end
    end
end

"""
    animateCar(track::Track, result, car::Car; fps=30, speedup=1.0, follow_car=false, view_radius=10.0)

Animate the car driving along the track.
- `follow_car=false` (default): global frame, full track visible.
- `follow_car=true`: camera follows the car, showing a window of `view_radius` around it.
"""
function animateCar(track::Track, result, car::Car; fps=30, speedup=1.0, follow_car=false, view_radius=10.0)
    xc, yc, thc, xl, yl, xr, yr = _precompute_track(track)

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())
    panel = _setup_panel(ax, track, car, follow_car, view_radius, xc, yc, thc, xl, yl, xr, yr)

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
function animateCarDual(track::Track, result, car::Car; fps=30, speedup=1.0, view_radius=10.0)
    xc, yc, thc, xl, yl, xr, yr = _precompute_track(track)

    fig = Figure(size=(1400, 700))
    ax_global = Axis(fig[1, 1], aspect=DataAspect(), title="Global")
    ax_follow = Axis(fig[1, 2], aspect=DataAspect(), title="Car Frame")

    panel_global = _setup_panel(ax_global, track, car, false, view_radius, xc, yc, thc, xl, yl, xr, yr)
    panel_follow = _setup_panel(ax_follow, track, car, true,  view_radius, xc, yc, thc, xl, yl, xr, yr)

    time_obs = Observable("t = 0.000 s")
    Label(fig[2, 1:2], time_obs; fontsize=20, halign=:center)
    rowsize!(fig.layout, 2, Fixed(30))

    display(GLMakie.Screen(), fig)

    _animate_loop!(fig, [panel_global, panel_follow], track, result, car, speedup, time_obs, xc, yc, xl, yl, xr, yr, view_radius)
    return fig
end
