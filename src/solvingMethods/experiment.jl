# ============================================================================
# Experiment configuration: disciplines, global constraints, solver backends,
# analysis toggles, and the top-level `Experiment` struct.
#
# Usage sketch:
#   exp = Experiment(
#       car = createTwintrack(true, track),
#       track = track,
#       discipline = Open(v_start=5.0),                 # or Closed()
#       solver = IpoptBackend(linear_solver="ma97"),    # or MadNLPBackend()
#       analysis = AnalysisConfig(animate=true, plot_jacobian=true),
#       global_constraints = [EnergyBudget(1.8e7)],     # Joules (optional)
#   )
#   run_experiment!(exp, segments, pol_order)
# ============================================================================


# ---------------------------------------------------------------------------
# Disciplines — define boundary conditions for the trajectory optimization.
# ---------------------------------------------------------------------------
abstract type Discipline end

"""
    Open(; v_start=5.0)

Open discipline: the car starts at a fixed longitudinal speed `v_start` (m/s)
and drives to the end of the track. Used for autocross, acceleration runs,
and any single-pass scenario.
"""
struct Open <: Discipline
    v_start::Float64
end
Open(; v_start::Real=5.0) = Open(float(v_start))

"""
    Closed()

Closed (periodic) discipline: the states at the start of the lap are equal to
the states at the end of the lap (except for the time coordinate which always
starts at zero). Used for endurance/skidpad-style closed-loop scenarios.
"""
struct Closed <: Discipline end


# ---------------------------------------------------------------------------
# Global constraints — lap-integral / global-level constraints on the model.
# Add new ones by subtyping GlobalConstraint and defining `apply_global!`.
# ---------------------------------------------------------------------------
abstract type GlobalConstraint end

"""
    EnergyBudget(max_energy_joules)

Limits the total drivetrain energy consumed over the run. Energy is
approximated as ∑ over collocation points of (Σ motor_torque × ω_wheel) × Δt,
using a trapezoidal rule on the path-to-time mapping. Only positive
(propulsive) torque contributes; regen is ignored in this simple version.

Use as an example of how to hook global constraints into the experiment.
"""
struct EnergyBudget <: GlobalConstraint
    max_energy::Float64  # Joules
end


# ---------------------------------------------------------------------------
# Solver backends — abstract over Ipopt / MadNLP. Each backend knows how to
# build a JuMP model configured for that solver.
# ---------------------------------------------------------------------------
abstract type SolverBackend end

"""
    IpoptBackend(; performSensitivity=false, attributes=Dict())

Ipopt configuration. Every Ipopt option lives in `attributes` so you can add or
remove keys without touching this struct. `performSensitivity=true` wraps the
model with DiffOpt for reverse-mode parameter sensitivities.

Example:
    IpoptBackend(attributes = Dict(
        "linear_solver"          => "ma97",
        "hsllib"                 => HSL_jll.libhsl_path,
        "max_iter"               => 3000,
        "tol"                    => 1e-6,
        "mu_strategy"            => "adaptive",
        "print_timing_statistics"=> "yes",
    ))
"""
Base.@kwdef mutable struct IpoptBackend <: SolverBackend
    performSensitivity::Bool = false
    attributes::Dict{String,Any} = Dict{String,Any}()
end

"""
    MadNLPBackend(; attributes=Dict())

MadNLP configuration. Every MadNLP option lives in `attributes`. For GPU:
    MadNLPBackend(attributes = Dict(
        "array_type"    => CUDA.CuArray,
        "linear_solver" => MadNLPGPU.CUDSSSolver,
    ))
"""
Base.@kwdef mutable struct MadNLPBackend <: SolverBackend
    attributes::Dict{String,Any} = Dict{String,Any}()
end


# ---------------------------------------------------------------------------
# Mesh refinement configuration.
# ---------------------------------------------------------------------------
"""
    MeshRefinementConfig(; tol, method, error_method, max_steps)

Controls the adaptive mesh refinement loop.

- `tol`:          segments with integrated error above this are split (default 1e-1)
- `method`:       `:h`  — insert midpoint into bad segments (refine mesh spacing)
                  `:p`  — increase polynomial order of bad segments (not yet impl.)
                  `:hp` — try p-refinement first, fall back to h-refinement
- `error_method`: passed to `getErrors`; `:ode` uses ODE residual, `:defect` uses collocation defect
- `max_iterations`: maximum number of refinement iterations (default 10)
"""
Base.@kwdef struct MeshRefinementConfig
    tol::Float64          = 1e-1
    method::Symbol        = :h        # :h | :p | :hp
    error_method::Symbol  = :ode      # :ode | :defect
    max_iterations::Int   = 10
end


# ---------------------------------------------------------------------------
# Analysis / post-processing toggles.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct AnalysisConfig
    plot_path::Bool = true
    plot_states::Bool = false
    plot_controls::Bool = false
    plot_jacobian::Bool = false
    plot_hessian::Bool = false
    animate::Bool = false
    animation_path::String = "results/animation.mp4"
    animation_speedup::Float64 = 1.0
    time_simulation::Bool = false                    # feedforward ODE verification
    plot_initialization::Bool = false
    save_results::Bool = false
    results_path::String = "results/"
end


# ---------------------------------------------------------------------------
# Experiment — the top-level container.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct Experiment
    car::Union{Nothing,Car} = nothing
    track::Union{Nothing,Track} = nothing
    discipline::Discipline = Open()
    solver::SolverBackend = IpoptBackend()
    analysis::AnalysisConfig = AnalysisConfig()
    global_constraints::Vector{GlobalConstraint} = GlobalConstraint[]
    mesh_refinement::MeshRefinementConfig = MeshRefinementConfig()
    # runtime state (populated by run_experiment!)
    model = nothing
    optiResult = nothing        # Result_interpolation, populated after solve
    params = nothing            # sensitivity parameter dict
end


# ---------------------------------------------------------------------------
# Model construction per backend.
# ---------------------------------------------------------------------------
function build_model(b::IpoptBackend)
    model = b.performSensitivity ?
        DiffOpt.nonlinear_diff_model(Ipopt.Optimizer) :
        Model(Ipopt.Optimizer)
    for (k, v) in b.attributes
        set_attribute(model, k, v)
    end
    return model
end

function build_model(b::MadNLPBackend)
    model = Model(MadNLP.Optimizer)
    for (k, v) in b.attributes
        set_attribute(model, k, v)
    end
    return model
end

# Convenience: does this backend want DiffOpt sensitivity?
wants_sensitivity(b::IpoptBackend) = b.performSensitivity
wants_sensitivity(::SolverBackend) = false


# ---------------------------------------------------------------------------
# Discipline boundary conditions.
# State ordering is assumed to be [vx, vy, ψ, ω, n, t] (6 states).
# ---------------------------------------------------------------------------
"""
    apply_boundary_conditions!(discipline, model, X, track)

Add the initial (and terminal, where applicable) state constraints to `model`.
"""
function apply_boundary_conditions!(d::Open, model::JuMP.Model, X, track)
    @constraint(model, X[1, 1] == d.v_start)           # vx
    @constraint(model, X[1, 2] == 0)                   # vy
    @constraint(model, X[1, 3] == track.theta[1])      # heading
    @constraint(model, X[1, 4] == 0)                   # ω
    @constraint(model, X[1, 6] == 0)                   # t
    return nothing
end

function apply_boundary_conditions!(::Closed, model::JuMP.Model, X, track)
    # Periodic: the car exits the lap in the same dynamic state it entered.
    # Time is not periodic — it's the objective.
    @constraint(model, X[1, 1] == X[end, 1])                  # vx
    @constraint(model, X[1, 2] == X[end, 2])                  # vy
    # Heading wrap: if the track is a closed loop, theta typically advances by 2π.
    Δθ = track.theta[end] - track.theta[1]
    @constraint(model, X[1, 3] == X[end, 3] - Δθ)             # ψ
    @constraint(model, X[1, 4] == X[end, 4])                  # ω
    @constraint(model, X[1, 5] == X[end, 5])                  # n (lateral)
    @constraint(model, X[1, 6] == 0)                          # t starts at 0
    return nothing
end


# ---------------------------------------------------------------------------
# Global constraint application.
# ---------------------------------------------------------------------------
"""
    apply_global!(constraint, model, X, U, car, s_all)

Add a lap-integral / global constraint to the model. Add new global
constraints by defining a new `apply_global!` method.
"""
function apply_global!(c::EnergyBudget, model::JuMP.Model, X, U, car::Car, s_all::AbstractVector)
    # Approximate energy: E ≈ Σ_i ΔT_i · P_i where P_i = Σ_motor τ_m · ω_wheel
    # ω_wheel = vx / radius (approx; ignores slip and vy)
    # ΔT_i = X[i+1, 6] - X[i, 6]
    # Motor torques live in U — we look them up via the control descriptor.
    #
    # This is a reasonable starter. For accurate energy budgeting you'd want
    # to include each motor/gearbox efficiency and the actual wheel angular
    # velocity, not vx/r.
    desc = car.carParameters.control_descriptor
    motor_control_idxs = Int[]
    for (idx, e) in enumerate(desc)
        if occursin("motor", lowercase(e.name)) && occursin("torque", lowercase(e.name))
            push!(motor_control_idxs, idx)
        end
    end
    if isempty(motor_control_idxs)
        @warn "EnergyBudget: no motor torque controls found in control_descriptor; skipping."
        return nothing
    end

    # Assume all drive tires share the same radius (reasonable for symmetric cars).
    r = car.drivetrain.tires[1].radius.value

    # Trapezoidal rule in time: E ≈ Σ_i 0.5 · (P_i + P_{i+1}) · Δt_i
    N = size(X, 1)
    # Per-node power expression.
    P_node = [sum(U[i, j] * X[i, 1] / r for j in motor_control_idxs) for i in 1:N]
    Δt = [X[i+1, 6] - X[i, 6] for i in 1:N-1]
    energy = sum(0.5 * (P_node[i] + P_node[i+1]) * Δt[i] for i in 1:N-1)

    @constraint(model, energy <= c.max_energy)
    return nothing
end


# ---------------------------------------------------------------------------
# Top-level driver.
# ---------------------------------------------------------------------------
"""
    run_experiment!(exp, segments, pol_order; variant="Lobatto")

Build the JuMP model, solve the trajectory optimization, and run whatever
post-processing is enabled in `exp.analysis`. Results land in `exp.result`.
"""
function run_experiment!(exp::Experiment, segments::Int64, pol_order::Int64; variant::String="Lobatto")
    exp.model = build_model(exp.solver)
    t_solve = @elapsed begin
        _, interp = find_optimal_trajectory_adaptive(exp, segments, pol_order, variant)
    end
    println("Solve time: $(round(t_solve, digits=2)) s")
    exp.optiResult = interp
    run_analysis!(exp)
    return exp
end

"""
    run_analysis!(exp)

Dispatch on `exp.analysis` flags to produce plots / animations / sensitivity.
"""
function run_analysis!(exp::Experiment)
    an = exp.analysis
    car   = exp.car
    track = exp.track
    res   = exp.optiResult

    if an.plot_path
        fig = Figure()
        ax = Axis(fig[1, 1], aspect=DataAspect())
        plotCarPath_interpolated(track, res, ax)
        display(GLMakie.Screen(), fig)

        if an.time_simulation
            sol = timeSimulation_interpolated(car, res, track)
            lines!(ax, getindex.(sol.u, 5), getindex.(sol.u, 6), label="Simulated in time")
            axislegend(ax, position=:rt)
        end
    end

    if an.plot_controls
        plot_controls_on_path(exp, res)
    end

    if an.plot_states
        plot_states_controls(car, res, track; sample_step=0.1)
    end

    if wants_sensitivity(exp.solver)
        sensitivityAnalysis(exp)
    end

    if an.animate
        animateCarDual(track, res, car;
            speedup = an.animation_speedup,
            view_radius = car.chassis.wheelbase.value * 3,
            cam_offset = 3.0,
            savepath = an.animation_path)
    end

    if an.plot_jacobian
        plot_jacobian_spy(exp.model)
    end
    if an.plot_hessian
        plot_hessian_spy(exp.model)
    end
    return nothing
end
