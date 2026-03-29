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
    animateCar(track::Track, result, car::Car; fps=30, speedup=1.0)

Animate the car driving along the track using simulation results.
Uses Observables for fast redraw — plot objects are created once, only data is updated.
"""
function animateCar(track::Track, result, car::Car; fps=30, speedup=1.0)
    fig = Figure(size=(1200, 800))
    ax = Axis(fig[1, 1], aspect=DataAspect(), title="Car Animation")
    plotTrack(track, b_plotStartEnd=false, ax=ax)

    # ── Create all Observables once ──
    wb = car.chassis.wheelbase.value
    tw = car.chassis.track.value
    dummy_pts = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)

    # Chassis
    chassis_obs = Observable(dummy_pts)
    poly!(ax, chassis_obs; color=(:steelblue, 0.3), strokecolor=:steelblue, strokewidth=2)

    # Tires (one Observable per wheel)
    tire_obs = [Observable(dummy_pts) for _ in car.wheelAssemblies]
    for obs in tire_obs
        poly!(ax, obs; color=:black, strokecolor=:gray30, strokewidth=1)
    end

    # Front wing
    front_wing_obs = Observable(dummy_pts)
    poly!(ax, front_wing_obs; color=:red, strokecolor=:darkred, strokewidth=1.5)

    # Rear wing
    rear_wing_obs = Observable(dummy_pts)
    poly!(ax, rear_wing_obs; color=:red, strokecolor=:darkred, strokewidth=1.5)

    # CoG marker
    cog_obs = Observable([Point2f(0, 0)])
    scatter!(ax, cog_obs; color=:yellow, markersize=8, marker=:cross, strokewidth=1.5, strokecolor=:black)

    # Trail
    trail_obs = Observable(Point2f[])
    lines!(ax, trail_obs; color=(:dodgerblue, 0.5), linewidth=2)

    display(GLMakie.Screen(), fig)

    # ── Precompute all frames ──
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

        ψ = state[3]
        n = state[5]
        steering = length(ctrl) >= 3 ? ctrl[3] : 0.0

        fc = track.fcurve(s)
        theta = fc[2]
        cx = fc[3] - n * sin(theta)
        cy = fc[4] + n * cos(theta)

        R = _rotmat2d(ψ)

        # Update chassis
        chassis_obs[] = _rect_points(cx, cy, wb + 0.2, tw + 0.1, ψ)

        # Update tires
        for (j, wa) in enumerate(car.wheelAssemblies)
            tire = car.drivetrain.tires[j]
            pos_local = wa.position.value[1:2]
            gp = R * pos_local .+ [cx, cy]
            steer = j <= 2 ? steering : 0.0
            tire_obs[j][] = _rect_points(gp[1], gp[2], 2 * tire.radius.value, tire.width.value * 0.8, ψ + steer)
        end

        # Update wings
        fc_pos = R * [wb / 2 + 0.15; 0.0] .+ [cx, cy]
        front_wing_obs[] = _rect_points(fc_pos[1], fc_pos[2], 0.05, tw + 0.3, ψ)
        rc_pos = R * [-wb / 2 - 0.15; 0.0] .+ [cx, cy]
        rear_wing_obs[] = _rect_points(rc_pos[1], rc_pos[2], 0.05, tw + 0.2, ψ)

        # Update CoG
        cog_obs[] = [Point2f(cx, cy)]

        # Update trail
        push!(trail_obs[], Point2f(cx, cy))
        notify(trail_obs)

        # Timing
        if i < n_frames
            t_current = state[6]
            if is_matrix
                t_next = result.states[i + 1, 6]
            else
                t_next = result.states(result.path[i + 1])[6]
            end
            sleep(max(0.0, (t_next - t_current) / speedup))
        end
    end

    return fig
end
