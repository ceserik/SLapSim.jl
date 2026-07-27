# =============================================================================
# compareVehicles.jl
#
# Renders a single video showing TWO solved vehicles racing on the same track at
# once — e.g. to compare different vehicle concepts (rear steering, torque
# vectoring, bigger/smaller accumulator, a different chassis altogether). Both
# cars are drawn with the full observable-based rendering (chassis outline,
# wheels + friction ellipses, aero), each tinted with its own colour halo + name
# label + trail so they stay visually distinguishable, and both are driven from
# a shared global clock so you can see who is ahead in real time.
#
# Usage
# -----
# Edit CAR_A / CAR_B below (factory functions of `track`, plus display names),
# then run:
#   julia --project=. src/experiments/compareVehicles.jl
#
# To compare arbitrary already-solved experiments from your own script, call
#   compare_cars(track, car1, res1, "Name 1", car2, res2, "Name 2"; output=...)
# directly instead of using the config block.
# =============================================================================

using Revise
using SLapSim
using CairoMakie
using Interpolations
using Printf

CairoMakie.activate!()

# ---------------------------------------------------------------------------
# Configuration — tweak these for the comparison you want.
# ---------------------------------------------------------------------------
const N_SEGMENTS = 30
const POL_ORDER  = 3
const FRAMERATE  = 30
const SPEEDUP    = 1.0                       # >1 = faster playback than real time
const OUTPUT     = "sync/animations/vehicle_comparison.mp4"
const HOLD_END   = 45                        # extra frames to linger once both cars finish

# The two vehicles to compare. Each factory takes `track` and returns a Car.
# Swap these for any two car-building functions to compare different concepts,
# e.g. `createBus`, or the same factory called with different design parameters.
CAR_A_NAME    = "Baseline"
CAR_A_FACTORY = track -> createTwintrack(true, track)
CAR_A_COLOR   = :dodgerblue

CAR_B_NAME    = "Variant"
CAR_B_FACTORY = track -> createTwintrack(true, track)
CAR_B_COLOR   = :crimson

track = doubleTurn(false, 0.1)

ipopt_attrs = Dict{String,Any}(
    "linear_solver" => "mumps",
    "alpha_for_y"   => "safer-min-dual-infeas",
    "print_level"   => 0,
)

function make_experiment(car)
    Experiment(
        car        = car,
        track      = track,
        discipline = Open(v_start = 5.0),
        solver     = IpoptBackend(performSensitivity = false, attributes = copy(ipopt_attrs)),
        mesh_refinement = MeshRefinementConfig(
            tol            = 1e-1,
            method         = :h,
            error_method   = :ode,
            max_iterations = 0,               # single solve on exactly N_SEGMENTS intervals
            segments       = N_SEGMENTS,
            pol_order      = POL_ORDER,
            variant        = "Radau",
        ),
        analysis = AnalysisConfig(
            plot_path = false, plot_states = false, plot_controls = false,
            plot_jacobian = false, plot_hessian = false, animate = false,
            time_simulation = false, plot_initialization = false,
        ),
    )
end

# ---------------------------------------------------------------------------
# Geometry helper: global pose (x, y, heading) of the car at arc-length s.
#   state[3] = heading psi,  state[5] = lateral offset n,  state[6] = time.
# ---------------------------------------------------------------------------
function _pose_at(res, track, s)
    st    = Float64.(res.states(s))
    fc    = track.fcurve(s)
    theta = fc[2]
    cx    = fc[3] - st[5] * sin(theta)
    cy    = fc[4] + st[5] * cos(theta)
    return Float64(cx), Float64(cy), Float64(st[3]), st
end

# Build a monotonic time -> arc-length map for a solved result, so a car's pose
# can be queried at any wall-clock time (held flat once it crosses the finish).
function _time_to_s(res, track)
    s_path = sort!(unique!(collect(Float64, res.path)))
    times  = [Float64(res.states(s)[6]) for s in s_path]
    times  = accumulate(max, times)                 # enforce monotonicity (time goes forward)
    keep   = [true; diff(times) .> 0]
    t_u, s_u = times[keep], s_path[keep]
    itp = extrapolate(interpolate((t_u,), s_u, Gridded(Linear())), Flat())
    return itp, t_u[end]                            # (time -> s), lap time T
end

function _free_path(path)
    isfile(path) || return path
    try
        rm(path; force = true)
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

# One persistent set of drawing observables + state for a single car on `ax`.
function _make_car_panel(ax, car, color, name)
    wb = car.chassis.wheelbase.value
    tw = car.chassis.track.value
    Fz_static = car.chassis.mass.value * 9.81 / length(car.wheelAssemblies)

    # Coloured halo behind the car so the two vehicles stay distinguishable
    # (chassis/aero outlines are transparent, tyres are always black).
    halo = Observable([Point2f(0, 0)])
    scatter!(ax, halo; color = (color, 0.25), markersize = 46, marker = :circle,
             strokewidth = 0)

    chassis_obs = car.chassis.setupObservables(ax)
    wa_obs      = [wa.setupObservables(ax, car.drivetrain.tires[i])
                   for (i, wa) in enumerate(car.wheelAssemblies)]
    aero_obs    = car.aero.setupObservables(ax)

    cog = Observable([Point2f(0, 0)])
    scatter!(ax, cog; color = color, marker = :cross, markersize = 10,
             strokewidth = 1.2, strokecolor = :black)

    trail = Observable(Point2f[])
    lines!(ax, trail; color = (color, 0.6), linewidth = 2)

    label_pos = Observable(Point2f(0, 0))
    text!(ax, label_pos; text = name, color = color, fontsize = 16,
          align = (:center, :bottom), offset = (0, 8))

    return (car = car, wb = wb, tw = tw, Fz_static = Fz_static,
            halo = halo, chassis = chassis_obs, wa = wa_obs, aero = aero_obs,
            cog = cog, trail = trail, label_pos = label_pos)
end

function _update_car_panel!(panel, track, res, s)
    car = panel.car
    cx, cy, psi, state = _pose_at(res, track, s)
    ctrl = Float64.(res.controls(s))

    car.stateMapping(state)
    car.controlMapping(ctrl)
    try
        car.carFunction(track, nothing)
    catch
        # Robust to any non-physical point (e.g. right at t=0); keep previous forces.
    end

    panel.halo[] = [Point2f(cx, cy)]
    car.chassis.updateObservables(panel.chassis, cx, cy, psi)
    for (j, wa) in enumerate(car.wheelAssemblies)
        wa.updateObservables(panel.wa[j], car.drivetrain.tires[j], cx, cy, psi, panel.Fz_static)
    end
    car.aero.updateObservables(panel.aero, cx, cy, psi, panel.wb, panel.tw)
    panel.cog[] = [Point2f(cx, cy)]
    panel.label_pos[] = Point2f(cx, cy)

    push!(panel.trail[], Point2f(cx, cy))
    notify(panel.trail)
    return nothing
end

# ---------------------------------------------------------------------------
# Public entry point: render two solved results racing side-by-side in time.
# ---------------------------------------------------------------------------
function compare_cars(track, car1, res1, name1::String, car2, res2, name2::String;
                      color1 = :dodgerblue, color2 = :crimson,
                      output = OUTPUT, framerate = FRAMERATE, speedup = SPEEDUP,
                      hold_end = HOLD_END)
    time2s_1, T1 = _time_to_s(res1, track)
    time2s_2, T2 = _time_to_s(res2, track)
    T_max = max(T1, T2)

    fig = Figure(size = (1920, 1080), backgroundcolor = :white)
    title_obs = Observable(_fmt_time(0.0))
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = :white,
              title = title_obs, titlesize = 22)

    plotTrack(track; ax = ax, b_plotStartEnd = false)

    panel1 = _make_car_panel(ax, car1, color1, name1)
    panel2 = _make_car_panel(ax, car2, color2, name2)

    legend_obs = Observable(_legend_text(name1, T1, name2, T2))
    Label(fig[2, 1], legend_obs; fontsize = 18, halign = :center)
    rowsize!(fig.layout, 2, Fixed(30))

    #resize_to_layout!(fig)

    n_frames = max(2, round(Int, framerate * (T_max / max(speedup, eps())))) + hold_end
    t_uniform = collect(range(0.0, T_max; length = n_frames - hold_end))

    output = _free_path(output)
    mkpath(dirname(output))
    println("Rendering comparison: $name1 (T=$(round(T1,digits=3))s) vs $name2 (T=$(round(T2,digits=3))s)")
    CairoMakie.record(fig, output, 1:n_frames; framerate = framerate) do fpos
        t = fpos <= length(t_uniform) ? t_uniform[fpos] : T_max
        s1 = time2s_1(t)
        s2 = time2s_2(t)
        _update_car_panel!(panel1, track, res1, s1)
        _update_car_panel!(panel2, track, res2, s2)
        title_obs[] = _fmt_time(t)
        legend_obs[] = _legend_text(name1, T1, name2, T2)
        print("\rRendering frame $fpos / $n_frames")
    end
    println("\nAnimation saved to $output")
    return fig
end

_fmt_time(t) = @sprintf("t = %.3f s", t)

function _legend_text(name1, T1, name2, T2)
    lead = T1 <= T2 ? name1 : name2
    gap  = abs(T1 - T2)
    return @sprintf("%s finishes at %.3f s   |   %s finishes at %.3f s   |   %s wins by %.3f s",
                     name1, T1, name2, T2, lead, gap)
end

# ---------------------------------------------------------------------------
# Solve both cars + render.
# ---------------------------------------------------------------------------
car1 = CAR_A_FACTORY(track)
car1.chassis.mass.value = 140.0
car2 = CAR_B_FACTORY(track)

exp1 = make_experiment(car1)
println("Solving $CAR_A_NAME …")
run_experiment!(exp1)

exp2 = make_experiment(car2)
println("Solving $CAR_B_NAME …")
run_experiment!(exp2)

fig = compare_cars(track, car1, exp1.optiResult, CAR_A_NAME, car2, exp2.optiResult, CAR_B_NAME;
                   color1 = CAR_A_COLOR, color2 = CAR_B_COLOR)

nothing
