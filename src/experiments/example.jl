using JuMP, Ipopt, Zygote
using Dierckx, DSP
using Interpolations
using RungeKutta
using  OrdinaryDiffEq
using DiffEqCallbacks
using DifferentiationInterface
using LinearAlgebra

# Type definitions (used by Motor, Tire, etc. structs)
const carVar = Union{Float64, JuMP.VariableRef, JuMP.AffExpr,
                     JuMP.NonlinearExpr, JuMP.QuadExpr}

mutable struct carParameter{T}
    value::T
    name::String
    unit::String
    size::Tuple{Int,Int}  # To store dimensions of the value (1,1) for scalar, (n,1) for vector, (n,m) for matrix
    role::Symbol  # :control, :state, :static parameter, :tunable parameter
    limits::Vector{Float64}
end

"""
I have created example that has the same problems as my bigger project. I used functions from my project and put them into single file the file, but it has 2000 lines :/ , I hope this form can be useful for you :sweat_smile:. 
The ODE function that describes the car dynamics is on line 1234. it is composed of different vehicle components, like motor, tire chassis, which define transformations of forces.

Collocation which creates constraints is on line 1505,



It really seems that the main issue is the number of terms inside constraints, I have 9 variables that describe the car model at each point at the track, the model contains  lot of transformations of forces between tires, wheel carries, CoG... So that may be too much :D
Should just add more variables to the model, or are there some other things that can increase the speed?

File is here, it works on my machine, pease let me know if it doesnt work for you.

https://github.com/ceserik/SLapSim.jl/blob/main/src/experiments/example.jl

these are all the packages required
```
using Pkg
Pkg.activate(".")
Pkg.add([
"JuMP",
"Ipopt",
"Zygote",
"Dierckx",
"DSP",
"Interpolations",
"RungeKutta",
"OrdinaryDiffEq",
"DiffEqCallbacks",
"DifferentiationInterface",
])
"""



# Component structs - defined early to avoid forward references
mutable struct Accumulator
    capacity::carParameter{carVar}
    maxPower::carParameter{carVar}
    minPower::carParameter{carVar}
    voltage::carParameter{carVar}
    current::carParameter{carVar}
    SoC::carParameter{carVar}
    mass::carParameter{carVar}
    resistance::carParameter{carVar}
end

mutable struct Gearbox{F,G}
    ratio::carParameter{carVar}
    torqueIn::carParameter{carVar}
    angularFrequencyIn::carParameter{carVar}
    torqueOut::carParameter{carVar}
    angularFrequencyOut::carParameter{carVar}
    loss::carParameter{carVar}
    compute::F
    setTorque::G
end

mutable struct Motor{F1,F2}
    torque::carParameter{carVar} #actual torque
    angularFrequency::carParameter{carVar} #actual angular freuency
    mass::carParameter{carVar}
    loss::carParameter{carVar} #actual power loss
    power::carParameter{carVar} # actual power draw
    torqueSpeedFunction::F1 #mapping speed to max torque
    constraints::F2
end

mutable struct Tire{F1,F2,F3,F4,F5}
    radius::carParameter{carVar}
    width::carParameter{carVar}
    inertia::carParameter{carVar}
    mass::carParameter{carVar}
    velocity::carParameter{Vector{carVar}}
    angularFrequency::carParameter{carVar}
    forces::carParameter{Vector{carVar}}
    slipAngle::carParameter{carVar}
    slipRatio::carParameter{carVar}
    compute::F1
    tireConstraints::F2
    setVelocity::F3
    maxSLipAngle::carParameter{carVar}
    scalingForce::carParameter{carVar}
    frictionCoefficient::carParameter{carVar}
    rollingResistance::carParameter{carVar}
    brakingForce::carParameter{carVar}
    setupObservables::F4
    updateObservables::F5
end

mutable struct Drivetrain
    motors::Vector{Motor}
    gearboxes::Vector{Gearbox}
    tires::Vector{Tire}
    accumulators::Accumulator
end
struct Result
    states::Matrix{Float64}
    controls::Matrix{Float64}
    path::Vector{Float64}
end

function _bary_eval(ref_nodes, weights, values, τ)
    for j in eachindex(ref_nodes)
        abs(τ - ref_nodes[j]) < 100 * eps(Float64) && return Float64(values[j])
    end
    num = sum(weights[j] * values[j] / (τ - ref_nodes[j]) for j in eachindex(ref_nodes))
    den = sum(weights[j] / (τ - ref_nodes[j]) for j in eachindex(ref_nodes))
    return num / den
end

function lessContraint(a,b,model=nothing)
    if isnothing(model)
        a = min(a,b)
    else
        @constraint(model,a<=b)
    end
    return a
end

function greaterContraint(a,b,model=nothing)
    if isnothing(model)
        a = max(a,b)
    else
        @constraint(model,a>=b)
    end
    return a
end

# Struct definitions needed by Car
struct VarEntry
    name::String
    targets::Vector{Pair}  # carParameter => vector_index (0 = scalar)
    lb::Float64
    ub::Float64
    role::Symbol  # :state or :control
end

mutable struct CarParameters
    mass::carParameter{carVar}
    inertia::carParameter{carVar}
    motorForce::carParameter{carVar}
    CL::carParameter{carVar}
    CD::carParameter{carVar}
    velocity::carParameter{Vector{carVar}}
    angularVelocity::carParameter{Vector{carVar}}
    psi::carParameter{carVar}
    n::carParameter{carVar}
    powerLimit::carParameter{carVar}
    lateralForce::carParameter{carVar}
    lateralTransfer::carParameter{carVar}
    brakeBias::carParameter{carVar}
    nControls::carParameter{carVar}
    nStates::carParameter{carVar}
    s::carParameter{carVar}
    state_descriptor::Vector{VarEntry}
    control_descriptor::Vector{VarEntry}
end

mutable struct Aero{F1,F2,F3}
    CL::carParameter{carVar}
    CD::carParameter{carVar}
    CoP::carParameter{carVar}
    compute::F1
    setupObservables::F2
    updateObservables::F3
end

mutable struct Chassis{F1,F2,F3}
    mass::carParameter{carVar}
    hitbox::F1
    wheelbase::carParameter{carVar}
    track::carParameter{carVar}
    CoG_X_pos::carParameter{carVar}
    CoG_Y_pos::carParameter{carVar}
    CoG_Z_pos::carParameter{carVar}
    width::carParameter{carVar}
    setupObservables::F2
    updateObservables::F3
end

mutable struct Suspension{F1,F2}
    tlong::carParameter{carVar}
    tlat::carParameter{carVar}
    theave::carParameter{carVar}
    stiffnessLong::carParameter{carVar}
    dampingLong::carParameter{carVar}
    stiffnessLat::carParameter{carVar}
    dampingLat::carParameter{carVar}
    stiffnessHeave::carParameter{carVar}
    dampingHeave::carParameter{carVar}
    lateralTransfer::carParameter{carVar}
    calculate::F1
    setInput ::F2
end

mutable struct WheelAssembly{F1,F2,F3,F4,F5}
    position::carParameter{Vector{carVar}}
    velocityPivot::carParameter{Vector{carVar}}
    velocityTire::carParameter{Vector{carVar}}
    maxAngle::carParameter{carVar}
    steeringAngle::carParameter{carVar}
    forces::carParameter{Vector{carVar}}
    torque::carParameter{Vector{carVar}}
    constraints::F1
    setVelocity::F2
    getTorque::F3
    setupObservables::F4
    updateObservables::F5
end

mutable struct Car
    carFunction::Function
    carParameters::CarParameters
    controlMapping::Function
    stateMapping::Function
    mapping::Function
    drivetrain::Drivetrain
    aero::Aero
    suspension::Suspension
    chassis::Chassis
    wheelAssemblies::Vector{WheelAssembly}
    state_desc::Union{Vector{VarEntry}, Nothing}
    control_desc::Union{Vector{VarEntry}, Nothing}
end

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
    sensitivity::Bool = false                        # run sensitivityAnalysis()
    save_results::Bool = false
    results_path::String = "results/"
end


# Track structs needed by Experiment
mutable struct Track
    curvature::Vector{Float64}
    rho::Vector{Float64}
    μ::Vector{Float64}
    sampleDistances::Vector{Float64}
    mapping::Function
    x::Vector{Float64}
    y::Vector{Float64}
    theta::Vector{Float64}
    widthR::Vector{Float64}
    widthL::Vector{Float64}
    inclination::Vector{Float64}
    slope::Vector{Float64}
    fcurve::Function
    s::Vector{Float64}
end

struct Track_interpolated
    s_nodes::Vector{Float64}
    x::Any
    y::Any
    heading::Any
    curvature::Any
    widthL::Any
    widthR::Any
    air_density::Any
    friction_c::Any
    inclination::Any
    fcurve::Any
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
        _, interp = find_optimal_trajectory_adaptive(exp, segments, pol_order, variant);
    end
    println("Solve time: $(round(t_solve, digits=2)) s");
    exp.optiResult = interp;
    run_analysis!(exp);
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
        #fig = Figure()
        #ax = Axis(fig[1, 1], aspect=DataAspect())
        #plotCarPath_interpolated(track, res, ax)
        #display(GLMakie.Screen(), fig)

        if an.time_simulation
            sol = timeSimulation_interpolated(car, res, track)
            #lines!(ax, getindex.(sol.u, 5), getindex.(sol.u, 6), label="Simulated in time")
            axislegend(ax, position=:rt)
        end
    end

    if an.plot_controls
#        plot_controls_on_path(exp, res)
    end

    if an.plot_states
     #   plot_states_controls(car, res, track; sample_step=0.1)
    end

    if an.sensitivity && wants_sensitivity(exp.solver)
        sensitivityAnalysis(exp)
    elseif an.sensitivity
        @warn "AnalysisConfig.sensitivity=true but solver backend does not support sensitivity; skipping."
    end

    if an.animate
        animateCarDual(track, res, car;
            speedup = an.animation_speedup,
            view_radius = car.chassis.wheelbase.value * 3,
            cam_offset = 3.0,
            savepath = an.animation_path)
    end

    if an.plot_jacobian
#        plot_jacobian_spy(exp.model)
    end
    if an.plot_hessian
 #       plot_hessian_spy(exp.model)
    end
    return nothing
end


function carParameter{T}(value::T, name::String, unit::String, role::Symbol=:static,limits::Vector{Float64} = [9999.0,-9999.0]) where T
    if isa(value, AbstractArray)
        size_tuple = (length(value), 1)  # For vectors, treat as (n,1)
        if ndims(value) == 2
            size_tuple = size(value)  # For matrices, use actual dimensions
        end
    else
        size_tuple = (1, 1)  # For scalars
    end

    return carParameter{T}(value, name, unit, size_tuple, role,limits)
end

# Convenience constructor for vector carParameters
carParameter{Vector{T}}(v::Vector{<:Real}, name::String, unit::String, args...) where T =
    carParameter{Vector{T}}(Vector{T}(v), name, unit, args...)

# Pretty printing for carParameter
function Base.show(io::IO, p::carParameter)
    val = p.value
    if isa(val, AbstractArray)
        print(io, "$(p.name) = $(val) [$(p.unit)] ($(p.role))")
    else
        print(io, "$(p.name) = $(val) [$(p.unit)] ($(p.role))")
    end
end

function Base.show(io::IO, ::MIME"text/plain", p::carParameter)
    println(io, "carParameter: $(p.name)")
    println(io, "  Value:  $(p.value)")
    println(io, "  Unit:   $(p.unit)")
    println(io, "  Size:   $(p.size)")
    println(io, "  Role:   $(p.role)")
    print(io,   "  Limits: $(p.limits)")
end

# Pretty printing for CarParameters
function Base.show(io::IO, ::MIME"text/plain", car::CarParameters)
    println(io, "CarParameters:")
    for fname in fieldnames(CarParameters)
        p = getfield(car, fname)
        println(io, "  $(p.name) = $(p.value) [$(p.unit)] ($(p.role))")
    end
end

# --- Variable descriptors for optimizer mapping ---


# Auto-read limits from the first target's carParameter
function VarEntry(name::String, targets::Vector{<:Pair}, role::Symbol=:state)
    param = first(targets).first
    VarEntry(name, collect(Pair, targets), param.limits[1], param.limits[2], role)
end

# Explicit bounds with typed pair vector
function VarEntry(name::String, targets::Vector{<:Pair}, lb::Real, ub::Real, role::Symbol=:state)
    VarEntry(name, collect(Pair, targets), Float64(lb), Float64(ub), role)
end

function apply_mapping!(desc::Vector{VarEntry}, values::AbstractVector)
    for (i, entry) in enumerate(desc)
        for (param, idx) in entry.targets
            if idx == 0
                param.value = values[i]
            else
                param.value[idx] = values[i]
            end
        end
    end
end

function apply_bounds!(desc::Vector{VarEntry})
    for entry in desc
        for (param, idx) in entry.targets
            val = idx == 0 ? param.value : param.value[idx]
            val isa Float64 || continue
            clamped = clamp(val, entry.lb, entry.ub)
            if idx == 0
                param.value = clamped
            else
                param.value[idx] = clamped
            end
        end
    end
end

get_bounds(desc::Vector{VarEntry}) = ([e.lb for e in desc], [e.ub for e in desc])

function get_scales(desc::Vector{VarEntry})
    return [max(abs(e.lb), abs(e.ub), 1.0) for e in desc]
end


Base.show(io::IO, ::MIME"text/plain", obj::Accumulator) = prettyPrintComponent(io, obj)


function createPepikCTU25()
    capacity = carParameter{carVar}(7.16,"ACP capacity","kWh")
    maxPower = carParameter{carVar}(80.0,"max Power out","kW")
    minPower = carParameter{carVar}(80.0,"min Power out","kW")
    voltage  = carParameter{carVar}(600.0,"Voltage","kW")
    current  = carParameter{carVar}(0.0,"Current","kW")
    SoC = carParameter{carVar}(100.0,"State of Charge","%")
    mass = carParameter{carVar}(38.49,"mass","kg")
    resistance = carParameter{carVar}(0.01,"Internal Resistance???","ohm")

    ACP = Accumulator(
        capacity,
        maxPower,
        minPower,
        voltage,
        current,
        SoC,
        mass,
        resistance
    )
    return ACP
end

#Base.show(io::IO, ::MIME"text/plain", obj) = prettyPrintComponent(io, obj)

const RHO_SEA_LEVEL = 1.225

function createBasicAero(; CL_a::Float64 = -5.0, CD_a::Float64 = 2.0, CoP_a::Float64 = 0.5)
    CL = carParameter{carVar}(CL_a,"Lift coeffcient","-")
    CD = carParameter{carVar}(CD_a,"Drag coeffcient","-")
    CoP = carParameter{carVar}(CoP_a,"Cener of pressure on front","-")
    function compute(vx::carVar, rho::Float64=RHO_SEA_LEVEL)
        downforce = -0.5 * rho * CL.value * vx^2
        drag = -0.5 * rho * CD.value * vx^2
        return (downforce=downforce, drag=drag)
    end
    function setupObservables(ax)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        front_obs = Observable(dummy)
        rear_obs  = Observable(dummy)
        poly!(ax, front_obs; color=:transparent, strokecolor=:transparent, strokewidth=1.5)
        poly!(ax, rear_obs;  color=:transparent, strokecolor=:transparent, strokewidth=1.5)
        return (front=front_obs, rear=rear_obs)
    end

    function updateObservables(obs::NamedTuple, x::Float64, y::Float64, ψ::Float64, wheelbase::Float64, trackwidth::Float64)
        R = _rotmat2d(ψ)
        fc = R * [wheelbase / 2 + 0.15; 0.0] .+ [x, y]
        obs.front[] = _rect_points(fc[1], fc[2], 0.05, trackwidth + 0.3, ψ)
        rc = R * [-wheelbase / 2 - 0.15; 0.0] .+ [x, y]
        obs.rear[] = _rect_points(rc[1], rc[2], 0.05, trackwidth + 0.2, ψ)
    end

    aero = Aero(
        CL,
        CD,
        CoP,
        compute,
        setupObservables,
        updateObservables,
    )
    return aero
end

function createCTU25chassis(; mass_p::Float64 = 280.0, wheelbase_p::Float64 = 1.525, track_p::Float64 = 1.2, CoGx::Float64 = 0.5, CogY::Float64 = 0.5, width_p::Float64 = 1.525)
    mass = carParameter{carVar}(mass_p,"mass","kg")
    wheelbase = carParameter{carVar}(wheelbase_p,"wheelbase","m")
    track = carParameter{carVar}(track_p,"track","m")
    width = carParameter{carVar}(width_p,"track","m")
    CoG_X_pos = carParameter{carVar}(CoGx,"ratio of CoG on front","-",:tunable)
    CoG_Y_pos = carParameter{carVar}(CogY,"ratio of CoG on left","-",:tunable)
    CoG_Z_pos = carParameter{carVar}(0.2,"CoG height in z from road","m",:tunable)
    function hitbox(n::carVar,racetrack::Union{Track,Nothing},model=nothing)
        #greaterContraint(n, -racetrack.widthR[1] + track.value / 2, model)
        #lessContraint(n, racetrack.widthL[1] - track.value / 2, model)
    end

    function setupObservables(ax)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        obs = Observable(dummy)
        poly!(ax, obs; color=:transparent, strokecolor=:transparent, strokewidth=0.5)
        return obs
    end

    function updateObservables(obs, x::Float64, y::Float64, ψ::Float64)
        obs[] = _rect_points(x, y, wheelbase.value + 0.2, track.value + 0.1, ψ)
    end

    chassis = Chassis(
        mass,
        hitbox,
        wheelbase,
        track,
        CoG_X_pos,
        CoG_Y_pos,
        CoG_Z_pos,
        width,
        setupObservables,
        updateObservables,
    )
    return chassis
end


Base.show(io::IO, ::MIME"text/plain", obj::Gearbox) = prettyPrintComponent(io, obj)

function createCTU25gearbox(ratio_p::Float64 = 11.46)
    ratio = carParameter{carVar}(ratio_p,"Gearbox ratio","-")
    torqueIn = carParameter{carVar}(0.0,"torque in","-")
    speedIn = carParameter{carVar}(0.0,"velocity in","rad/s")

    torqueOut = carParameter{carVar}(0.0,"torque out","-")
    angularFrequencyOut = carParameter{carVar}(0.0,"velocity out","rad/s")
    loss = carParameter{carVar}(0.0,"loss","Watt")

    function compute()
        torqueOut.value = torqueIn.value * ratio.value
        angularFrequencyOut.value = speedIn.value / ratio.value
    end

    function setTorque(torque::carVar)
        torqueIn.value = torque
    end
    
    gearbox = Gearbox(
        ratio,
        torqueIn,
        speedIn,
        torqueOut,
        angularFrequencyOut,
        loss,
        compute,
        setTorque
    )
end


Base.show(io::IO, ::MIME"text/plain", obj::Motor) = prettyPrintComponent(io, obj)

function createFischerMotor(maxTorqueVal::Float64=29.0,maxTorqueVal_regen::Float64=29.0)
    torque = carParameter{carVar}(0.0,"motor torque","Nm",:static,[-maxTorqueVal_regen,maxTorqueVal])
    angularFrequency = carParameter{carVar}(0.0,"angular frequency","rad/s")
    loss = carParameter{carVar}(0.0,"loss","W")
    torqueSpeedFunction = angularFrequency::Float64 -> maxTorqueVal
    mass = carParameter{carVar}(3.0,"mass??","kg")
    power = carParameter{carVar}(3.0,"Power draw","kW")

    function constraints(u,model=nothing)
        maxTorque = torqueSpeedFunction(0.0)
#        u = lessContraint(u/maxTorque, 1.0, model) * maxTorque
#        u = greaterContraint(u/maxTorque, -1.0, model) * maxTorque
        return u
    end


    motor = Motor(
        torque,
        angularFrequency,
        mass,
        loss,
        power,
        torqueSpeedFunction,
        constraints
    )
    return motor
end

Base.show(io::IO, ::MIME"text/plain", obj::Suspension) = prettyPrintComponent(io, obj)

function createSimpleSuspension()
    tlong = carParameter{carVar}(0.0,"longitudinal transfer time constant","s")
    tlat = carParameter{carVar}(0.0,"lateral transfer time constant","s")
    theave = carParameter{carVar}(0.0,"heave transfer time constant","s")

    stiffnessLong = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    stiffnessLat = carParameter{carVar}(0.0,"Lat stifness","N/rad??")
    stiffnessHeave = carParameter{carVar}(0.0,"Long stifness","N/rad??")
    
    dampingLong = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingLat = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    dampingHeave = carParameter{carVar}(0.0,"long dampnng","N/rad/s`")
    lateralTransfer = carParameter{carVar}(0.0,"lateral load transfer","N")
    chassis::Union{Chassis, Nothing} = nothing


    function setInput(chassis_in::Chassis)
        chassis = chassis_in

    end

    function calculate(downforce::carVar=0.0, CoP::carVar=0.5)
        totalLoad = chassis.mass.value * 9.81 + downforce
        frontRatio = chassis.CoG_X_pos.value
        rearRatio = 1 - chassis.CoG_X_pos.value
        leftRatio = 1 - chassis.CoG_Y_pos.value
        rightRatio = chassis.CoG_Y_pos.value
        weightFront = chassis.mass.value * 9.81 * frontRatio
        weightRear = chassis.mass.value * 9.81 * rearRatio
        aeroFront = downforce * CoP
        aeroRear = downforce * (1 - CoP)
        Fz_FL = (weightFront + aeroFront) * leftRatio
        Fz_FR = (weightFront + aeroFront) * rightRatio
        Fz_RL = (weightRear + aeroRear) * leftRatio
        Fz_RR = (weightRear + aeroRear) * rightRatio
        return [Fz_FL, Fz_FR, Fz_RL, Fz_RR]
    end

    susp = Suspension(
        tlong,
        tlat,
        theave,
        stiffnessLong,
        dampingLong,
        stiffnessLat,
        dampingLat,
        stiffnessHeave,
        dampingHeave,
        lateralTransfer,
        calculate,
        setInput    )
    return susp
end

const _N_ELLIPSE_PTS = 41
const _ELLIPSE_ANGLES = range(0, 2π, length=_N_ELLIPSE_PTS)

mutable struct  Tire{F1,F2,F3,F4,F5}
    radius::carParameter{carVar}
    width::carParameter{carVar}
    inertia::carParameter{carVar}
    mass::carParameter{carVar}
    velocity::carParameter{Vector{carVar}}
    angularFrequency::carParameter{carVar}
    forces::carParameter{Vector{carVar}}
    slipAngle::carParameter{carVar}
    slipRatio::carParameter{carVar}
    compute::F1
    tireConstraints::F2
    setVelocity::F3
    maxSLipAngle::carParameter{carVar}
    scalingForce::carParameter{carVar}
    frictionCoefficient::carParameter{carVar}
    rollingResistance::carParameter{carVar}
    brakingForce::carParameter{carVar}
    setupObservables::F4
    updateObservables::F5
end

function createR20_pacejka(motor::Motor,gearbox::Gearbox)
    radius = carParameter{carVar}(0.205, "tire radius", "m")
    width = carParameter{carVar}(0.3, "tire width, wrong", "m")
    inertia = carParameter{carVar}(0.3, "tire inertia, wrong", "m")
    mass = carParameter{carVar}(1.0, "tire mass,wrong", "kg")
    velocity = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "velocity", "m/s")
    angularFrequency = carParameter{carVar}(0.0, "angular velocity", "rad/s")
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Tire Force", "N")
    slipAngle = carParameter{carVar}(0.0, "Slip angle", "rad")
    slipRatio = carParameter{carVar}(0.0, "slip ratio", "-")
    maxSlipAngle = carParameter{carVar}(5/180*pi, "max slip angle", "rad")
    scalingForce = carParameter{carVar}(0.0, "max longitudinal force", "N")
    scalingForce.value = motor.torqueSpeedFunction(0.0) * gearbox.ratio.value / radius.value
    frictionCoefficient = carParameter{carVar}(1.0, "Friction Coefficient", "-", :tunable)
    rollingResistance = carParameter{carVar}(-0.010, "Rolling resistance coeff", "-", :tunable)
    brakingForce = carParameter{carVar}(0.0, "braking force", "N")
    # Pacejka Magic Formula coefficients
    B = 9.62   # stiffness factor
    C = 2.59  # shape factor
    E = 1  # curvature factor

    function setVelocity(velocityIn::Vector{carVar})
        velocity.value = velocityIn
    end

    function compute(inTorque::carVar, optiModel::Union{JuMP.Model,Nothing}=nothing)
        slipAngle.value = -atan(velocity.value[2], velocity.value[1])
        α = slipAngle.value * 1
        D = forces.value[3] * frictionCoefficient.value
        forces.value[2] = D * sin(C * atan(B * α - E * (B * α - atan(B * α))))
        forces.value[1] = inTorque/radius.value + brakingForce.value + forces.value[3] * rollingResistance.value
    end

    function tireConstraints(model=nothing)
        μFz = frictionCoefficient.value * forces.value[3]
        fx = forces.value[1] / μFz
        fy = forces.value[2] / μFz
        lessContraint(fx^2 + fy^2, 1.1, model) #crucial to remove oscilatory behavior with bus, should look deeper into tire model suitable for optimisation
    end

    function setupObservables(ax)
        dummy = _rect_points(0.0, 0.0, 1.0, 1.0, 0.0)
        rect_obs = Observable(dummy)
        poly!(ax, rect_obs; color=(:black, 0.3), strokecolor=(:gray30, 0.5), strokewidth=1)
        ellipse_obs = Observable([Point2f(0, 0) for _ in _ELLIPSE_ANGLES])
        lines!(ax, ellipse_obs; color=(:gray, 0.5), linewidth=1, linestyle=:dash)
        force_pos_obs = Observable([Point2f(0, 0)])
        force_dir_obs = Observable([Point2f(0, 0)])
        arrows2d!(ax, force_pos_obs, force_dir_obs; color=:green, shaftwidth=4, tipwidth=8, tiplength=8)
        return (rect=rect_obs, ellipse=ellipse_obs, force_pos=force_pos_obs, force_dir=force_dir_obs)
    end

    function updateObservables(obs::NamedTuple, wx::Float64, wy::Float64, θ::Float64, Fz_static::Float64)
        r = radius.value
        s = 4r / Fz_static
        R = _rotmat2d(θ)
        obs.rect[] = _rect_points(wx, wy, 2r, width.value * 0.8, θ)
        Fmax_s = abs(forces.value[3]) * frictionCoefficient.value * s
        obs.ellipse[] = [Point2f(R * [Fmax_s * cos(a), Fmax_s * sin(a)] + [wx, wy]) for a in _ELLIPSE_ANGLES]
        f_global = R * [forces.value[1], forces.value[2]] .* s
        obs.force_pos[] = [Point2f(wx, wy)]
        obs.force_dir[] = [Point2f(f_global[1], f_global[2])]
    end

    tire = Tire(
        radius,
        width,
        inertia,
        mass,
        velocity,
        angularFrequency,
        forces,
        slipAngle,
        slipRatio,
        compute,
        tireConstraints,
        setVelocity,
        maxSlipAngle,
        scalingForce,
        frictionCoefficient,
        rollingResistance,
        brakingForce,
        setupObservables,
        updateObservables,
    )
end

Base.show(io::IO, ::MIME"text/plain", obj::WheelAssembly) = prettyPrintComponent(io, obj)

function createBasicWheelAssembly(position::Vector{carVar})
    maxAngle = carParameter{carVar}(25 / 180 * pi, " max steering angle", "rad")
    steeringAngle = carParameter{carVar}(0.0, "steering angle", "rad", :control, [-maxAngle.value, maxAngle.value])
    forces = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Pivot forces", "N N N")
    torque = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Generated torque on CoG", "N N N")
    velocityPivot = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Velocity at pivot", "m/s")
    velocityTire = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Velocity in tire frame", "m/s")
    position = carParameter{Vector{carVar}}(position, "Position from CoG", "m m m")
    function rotZ(steering::carVar) #check direction
        out = [
            cos(steering) -sin(steering) 0;
            sin(steering) cos(steering) 0;
            0 0 1
        ]
        return out
    end

    function rotZinv(steering::carVar) #check direction
        out = [
            cos(steering) sin(steering) 0;
            -sin(steering) cos(steering) 0;
            0 0 1
        ]
        return out
    end
    function constraints(model=nothing)
        if !isnothing(model)
        #    set_lower_bound(steeringAngle.value, -maxAngle.value)
        #    set_upper_bound(steeringAngle.value,  maxAngle.value)
        end
    end

    function pivot2CoG(forces::Vector{carVar})
        #returns moments on cog from wheel forces
        torque.value = cross(position.value, forces)
    end

    function setPivotForce(forcesIn::Vector{carVar})
        forces.value = rotZ(steeringAngle.value) * forcesIn
    end

    function setPivotVelocity(angularVelocity::Vector{carVar}, CoGvelocity::Vector{carVar})
        velocityPivot.value = CoGvelocity + cross(angularVelocity, position.value) #CoG2wheelAssembly(velocity.value,angularVelocity)
    end

    function setTireSpeeds()
        velocityTire.value = rotZinv(steeringAngle.value) * velocityPivot.value
    end

    function setVelocity(angularVelocity::Vector{carVar}, velocity::Vector{carVar})
        setPivotVelocity(angularVelocity, velocity)
        setTireSpeeds()
    end

    function getTorque(forcesIn::Vector{carVar})
        setPivotForce(forcesIn)
        pivot2CoG(forces.value)
    end

    function setupObservables(ax, tire::Tire)
        return tire.setupObservables(ax)
    end

    function updateObservables(obs::NamedTuple, tire::Tire, x::Float64, y::Float64, ψ::Float64, Fz_static::Float64)
        R = _rotmat2d(ψ)
        global_pos = R * [position.value[1], position.value[2]] + [x, y]
        total_angle = ψ + steeringAngle.value
        tire.updateObservables(obs, global_pos[1], global_pos[2], total_angle, Fz_static)
    end

    testWheelAssembly = WheelAssembly(
        position,
        velocityPivot,
        velocityTire,
        maxAngle,
        steeringAngle,
        forces,
        torque,
        constraints,
        setVelocity,
        getTorque,
        setupObservables,
        updateObservables,
    )
    return testWheelAssembly
end

function createTwintrack(pacejka::Bool=true,track::Union{Track,Nothing} = nothing)

    if isnothing(track)
        widthR = 1.5
        widthL = 1.5
    else
        # Use full-track envelope so the static n limits are the loosest
        # legal bounds anywhere on the track. Adaptive solver will tighten
        # them per node from Track_interpolated. Margin kept at 0.6.
        widthL = maximum(track.widthL)
        widthR = maximum(track.widthR)
    end

    velocity         = carParameter{Vector{carVar}}([10.0, 10.0, 0.0], "Velocity", "m/s", :static, [3.0, 60.0])
    angularVelocity  = carParameter{Vector{carVar}}([0.0, 0.0, 0.0], "Angular Velocity", "rad/s", :static, [-10.0, 10.0])

    mass             = carParameter{carVar}(280.0, "Mass", "kg", :tunable, [200.0, 320.0])
    motorForce       = carParameter{carVar}(1000.0, "motorForce", "N")
    lateralForce     = carParameter{carVar}(0.0, "lateral Force", "N")
    lateralTransfer  = carParameter{carVar}(0.0, "lateral load transfer", "N")
    brakeBias        = carParameter{carVar}(0.6, "brake bias front", "-", :tunable)
    CL               = carParameter{carVar}(5.0, "Lift Coefficient", "-")
    CD               = carParameter{carVar}(2.0, "Drag Coefficient", "-")
    powerLimit       = carParameter{carVar}(80000.0, "PowerLimit", "W")
    psi              = carParameter{carVar}(0.0, "heading", "rad", :static, [-4pi, 4pi])
    n                = carParameter{carVar}(0.0, "Distance from centerline", "m", :static, [-widthR+0.6, widthL-0.6])
    brakeCommand     = carParameter{carVar}(0.0, "brake command", "N", :control, [-20000.0, 0.0])
    nControls        = carParameter{carVar}(0.0, "number of controlled parameters", "-")
    inertia          = carParameter{carVar}(100.0, "Inertia", "kg*m^2", :tunable)
    nStates          = carParameter{carVar}(0.0, "number of car states", "-")
    s                = carParameter{carVar}(1.0, "longitudinal position on track", "-", :static, [0.0, 200.0])

    gearboxFL = createCTU25gearbox()
    gearboxFR = createCTU25gearbox()
    gearboxRL = createCTU25gearbox()
    gearboxRR = createCTU25gearbox()

    motorFL = createFischerMotor()
    motorFR = createFischerMotor()
    motorRL = createFischerMotor()
    motorRR = createFischerMotor()

    if pacejka
        tireFL = createR20_pacejka(motorFL, gearboxFL)
        tireFR = createR20_pacejka(motorFR, gearboxFR)
        tireRL = createR20_pacejka(motorRL, gearboxRL)
        tireRR = createR20_pacejka(motorRR, gearboxRR)
    else
        tireFL = createR20lin(motorFL, gearboxFL)
        tireFR = createR20lin(motorFR, gearboxFR)
        tireRL = createR20lin(motorRL, gearboxRL)
        tireRR = createR20lin(motorRR, gearboxRR)

    end
    drivetrain = Drivetrain(
        [motorFL, motorFR, motorRL, motorRR],
        [gearboxFL, gearboxFR, gearboxRL, gearboxRR],
        [tireFL, tireFR, tireRL, tireRR],
        createPepikCTU25())

    aero = createBasicAero()
    suspension = createSimpleSuspension()

    chassis = createCTU25chassis()
    suspension.setInput(chassis)
    cogOffsetX = (chassis.CoG_X_pos.value - 0.5) * chassis.wheelbase.value
    cogOffsetY = (chassis.CoG_Y_pos.value - 0.5) * chassis.track.value
    wheelAssemblyFL = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyFR = createBasicWheelAssembly(Vector{carVar}([chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRL = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, chassis.track.value / 2 - cogOffsetY, 0]))
    wheelAssemblyRR = createBasicWheelAssembly(Vector{carVar}([-chassis.wheelbase.value / 2 - cogOffsetX, -chassis.track.value / 2 - cogOffsetY, 0]))

    

    wheelAssemblies = [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL, wheelAssemblyRR]

    state_descriptor = VarEntry[
        VarEntry("vx",    [velocity => 1],          :state),
        VarEntry("vy",    [velocity => 2], -1.5, 1.5, :state),
        VarEntry("psi",   [psi => 0],               :state),
        VarEntry("omega", [angularVelocity => 3],   :state),
        VarEntry("n",     [n => 0],                 :state),
        VarEntry("t",     [s => 0],                 :state),
    ]

    control_descriptor = VarEntry[
        VarEntry("torque_rear",  [drivetrain.motors[3].torque => 0, drivetrain.motors[4].torque => 0], :control),
        VarEntry("steering",    [wheelAssemblies[1].steeringAngle => 0, wheelAssemblies[2].steeringAngle => 0], :control),
        VarEntry("torque_front", [drivetrain.motors[1].torque => 0, drivetrain.motors[2].torque => 0], :control),
        VarEntry("brake", [brakeCommand => 0], :control),
    ]
#    @infiltrate
    nControls.value = Float64(length(control_descriptor))
    nStates.value   = Float64(length(state_descriptor))

    function controlMapping(controls::AbstractVector)
        apply_mapping!(control_descriptor, controls)
        apply_bounds!(control_descriptor)
        frontBrake = brakeCommand.value * brakeBias.value
        rearBrake = brakeCommand.value * (1 - brakeBias.value)
        drivetrain.tires[1].brakingForce.value = frontBrake
        drivetrain.tires[2].brakingForce.value = frontBrake
        drivetrain.tires[3].brakingForce.value = rearBrake
        drivetrain.tires[4].brakingForce.value = rearBrake
    end

    function stateMapping(states::AbstractVector)
        apply_mapping!(state_descriptor, states)
    end

    function anyTrack(track::Union{Track,Nothing}=nothing, optiModel::Union{JuMP.Model,Nothing}=nothing)
        # Transformation of velocities from cog to wheels and steering
        if !isnothing(optiModel)
            # Prevent simultaneous acceleration and braking
            # If torque > 0 (accelerating), brake must be 0; if brake < 0 (braking), torque must be 0
            T_max = drivetrain.motors[3].torque.limits[2]  # 29 Nm
            F_max = -brakeCommand.limits[1]                # 20000 N
            # Product complementarity: (torque/T_max) * (-brake/F_max) ≤ ε
            # Both factors are in [0,1], so product near zero enforces mutual exclusion
            @constraint(optiModel, (drivetrain.motors[3].torque.value / T_max) * (-brakeCommand.value / F_max) <= 0.001)
            @constraint(optiModel, (drivetrain.motors[1].torque.value / T_max) * (-brakeCommand.value / F_max) <= 0.001)
        end
        for wa in wheelAssemblies
            wa.setVelocity(angularVelocity.value, velocity.value)
            wa.constraints(nothing)
        end

        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].setVelocity(wheelAssemblies[i].velocityTire.value)
        end

        #gearing of forces from motor to tire
        for i in eachindex(drivetrain.gearboxes)
            drivetrain.gearboxes[i].setTorque(drivetrain.motors[i].torque.value)
            drivetrain.gearboxes[i].compute()
        end



        aeroForces = aero.compute(velocity.value[1], isnothing(track) ? RHO_SEA_LEVEL : track.rho[1])

        forces = suspension.calculate(aeroForces.downforce, aero.CoP.value)
        for i in eachindex(forces)
            drivetrain.tires[i].forces.value[3] = forces[i]
        end

        #tire forces
        for i in eachindex(drivetrain.tires)
            drivetrain.tires[i].compute(drivetrain.gearboxes[i].torqueOut.value, optiModel)
        end

        ######################################################CONSTRAINTS###########################
        # motor torque limit — one per shared control
        #drivetrain.motors[1].constraints(drivetrain.motors[1].torque.value, optiModel) # front (controls[3])
        #drivetrain.motors[3].constraints(drivetrain.motors[3].torque.value, optiModel) # rear  (controls[1])
        #steering angle
        #wheelAssemblies[1].constraints(optiModel)
        #hitbox
        #chassis.hitbox(n.value, track, optiModel)
        for tire in drivetrain.tires
            tire.tireConstraints(optiModel)
        end

        for i in eachindex(wheelAssemblies)

            wheelAssemblies[i].getTorque(drivetrain.tires[i].forces.value)
        end


        cogForce = zero(wheelAssemblies[1].forces.value)
        cogMoment = zero(wheelAssemblies[1].forces.value)
        for i in eachindex(wheelAssemblies)
            cogForce = cogForce .+ wheelAssemblies[i].forces.value
            cogMoment = cogMoment + wheelAssemblies[i].torque.value
        end
        cogForce = cogForce .+ [aeroForces.drag, 0.0, 0.0]

        dv = cogForce ./ mass.value - angularVelocity.value × velocity.value #really check what sign should be here !!!! podla mna bednarik skripta fyzika1 Kapitola 8  Neinerciální vztažné soustavy, neboli pro vyjádření časové změny libovolné vektorové veličiny v nečárkované soustavě je možné použít následujícího operátoru:  d·  dt = d′·  dt + ω × · , (8.11)
        dangularVelocity = cogMoment ./ inertia.value
        dx = [dv[1], dv[2], angularVelocity.value[3], dangularVelocity[3]]
        return dx

    end

    p = CarParameters(
        mass,
        inertia,
        motorForce,
        CL,
        CD,
        velocity,
        angularVelocity,
        psi,
        n,
        powerLimit,
        lateralForce,
        lateralTransfer,
        brakeBias,
        nControls,
        nStates,
        s,
        state_descriptor,
        control_descriptor
    )

    afto = Car(
        anyTrack,
        p,
        controlMapping,
        stateMapping,
        f -> (0.0),
        drivetrain,
        aero,
        suspension,
        chassis,
        [wheelAssemblyFL, wheelAssemblyFR, wheelAssemblyRL, wheelAssemblyRR],
        state_descriptor,
        control_descriptor,
    )
    return afto
end



function interp1(X, V, Xq,type=nothing)
    if type == "PCHIP"
        spl = Spline1D(X, V)
        return spl(Xq)
    else
        knots = (X,)
        itp = interpolate(knots, V, Gridded(Linear()))
        itp.(Xq)
    end
end
struct Result_interpolation
    states
    controls
    path
    segment_edges
end

struct Collocation
    createConstraints
    createInterpolator
    tableau
    interpolator::Union{Function,Nothing}
    f
end

function barycentric_weights(x)
    w = zeros(length(x))

    for j = eachindex(x)
        tmp = 1
        for k = eachindex(x)
            if k != j
                tmp *= (x[j] - x[k])
            end
        end
        w[j] = 1 / tmp
    end
    return w
end



function make_result_interpolation(x::AbstractMatrix{Float64}, u::AbstractMatrix{Float64}, s_states::Vector{Float64}, s_controls::Vector{Float64})
    state_interps = [linear_interpolation(s_states, x[:, i]) for i in 1:size(x, 2)]
    control_interps = [extrapolate(interpolate((s_controls,), u[:, i], Gridded(Linear())), Flat()) for i in 1:size(u, 2)]

    states_interp(t) = [itp(t) for itp in state_interps]
    controls_interp(t) = [itp(t) for itp in control_interps]
    return Result_interpolation(states_interp, controls_interp, s_states, Float64[])
end

make_result_interpolation(x::AbstractMatrix{Float64}, u::AbstractMatrix{Float64}, s::Vector{Float64}) = make_result_interpolation(x, u, s, s)

function trackMapping(track::Track,trackCopy::Track ,s)
    trackCopy.curvature = interp1(track.sampleDistances,track.curvature,s)
    trackCopy.theta = interp1(track.sampleDistances,track.theta,s)
end

function createLobattoIIIA_Adaptive(f, stages, model, nControls, nStates, track;
                                    x_scale=ones(nStates), u_scale=ones(nControls),
                                    x_lb=s -> -x_scale, x_ub=s -> x_scale,
                                    u_lb=s -> -u_scale, u_ub=s -> u_scale)
    tableau = TableauLobattoIIIA(stages)
    τ_ref  = tableau.c
    w_bary = barycentric_weights(τ_ref)

    # Evaluate a bound function at path position and verify that it returns the expected number of entries.
    function _eval_bound(bound_fn, s, expected_length)
        bound_values = bound_fn(s)
        length(bound_values) == expected_length || error(
            "bound function returned vector of length $(length(bound_values)), " *
            "expected $(expected_length) at s=$(s)")
        return bound_values
    end

    function createDynamicConstraints(segment_edges, initialization)
        number_of_segments = length(segment_edges) - 1
        totalPoints = number_of_segments * (stages - 1) + 1

        # Build s_all FIRST so per-node bounds can be evaluated during variable creation.
        s_all = zeros(totalPoints)
        h_all = diff(segment_edges)
        idx = 1
        for segment = 1:number_of_segments #boundary nodes
            s_all[idx] = segment_edges[segment]
            h = segment_edges[segment+1] - segment_edges[segment]
            idx += 1

            for point = 2:(stages-1) #inner nodes
                #@infiltrate
                s_all[idx] = segment_edges[segment] + h * tableau.c[point]
                idx += 1
            end
        end
        s_all[end] = segment_edges[end]

        X_raw = Matrix{VariableRef}(undef, totalPoints, nStates)
        U_raw = Matrix{VariableRef}(undef, totalPoints, nControls)

        # Print bounds & scales summary using the value at s=0 (representative).
        let s0 = s_all[1]
            x_lb0 = _eval_bound(x_lb, s0, nStates)
            x_ub0 = _eval_bound(x_ub, s0, nStates)
            u_lb0 = _eval_bound(u_lb, s0, nControls)
            u_ub0 = _eval_bound(u_ub, s0, nControls)
            println("State bounds & scales (at s=$(round(s0,digits=3))):")
            for j in eachindex(x_lb0)
                println("  x[$j]: lb=$(round(x_lb0[j],digits=3))  ub=$(round(x_ub0[j],digits=3))  scale=$(round(x_scale[j],digits=3))")
            end
            println("Control bounds & scales (at s=$(round(s0,digits=3))):")
            for j in eachindex(u_lb0)
                println("  u[$j]: lb=$(round(u_lb0[j],digits=3))  ub=$(round(u_ub0[j],digits=3))  scale=$(round(u_scale[j],digits=3))")
            end
        end

        for i = 1:totalPoints
            si = s_all[i]
            x_lb_i = _eval_bound(x_lb, si, nStates) ./ x_scale
            x_ub_i = _eval_bound(x_ub, si, nStates) ./ x_scale
            u_lb_i = _eval_bound(u_lb, si, nControls) ./ u_scale
            u_ub_i = _eval_bound(u_ub, si, nControls) ./ u_scale
            for j = 1:nStates
                X_raw[i, j] = @variable(model, lower_bound=x_lb_i[j], upper_bound=x_ub_i[j])
            end
            for j = 1:nControls
                U_raw[i, j] = @variable(model, lower_bound=u_lb_i[j], upper_bound=u_ub_i[j])
            end
        end

        #scaling
        X = X_raw .* x_scale'
        U = U_raw .* u_scale'

        # ============================ Odow test ==========================================
        #F_raw = Matrix{VariableRef}(undef, totalPoints, nStates)
        #for i = 1:totalPoints, j = 1:nStates
        #    F_raw[i, j] = @variable(model)
        #end
        #for i = 1:totalPoints
        #    dx_i = f(X[i, :], U[i, :], s_all[i], model)
        #    @constraint(model, F_raw[i, :] .== dx_i ./ x_scale)
        #end
        #segment_start_idx = 1
        #for segment = 1:number_of_segments
        #    for stage = 2:stages
        #        fxSum = sum(tableau.a[stage, col] .* F_raw[segment_start_idx+col-1, :] for col in 1:stages)
        #        @constraint(model, X_raw[segment_start_idx+stage-1, :] .== X_raw[segment_start_idx, :] .+ h_all[segment] .* fxSum)
        #    end
        #    segment_start_idx = segment_start_idx + stages - 1
        #end


        # === end intermediate-F variant ===

        # ==============================Original =================================


        segment_start_idx = 1
        f(X[segment_start_idx, :], U[segment_start_idx, :], s_all[segment_start_idx], model)
        for segment = 1:number_of_segments
            F = [f(X[segment_start_idx+col-1, :], U[segment_start_idx+col-1, :], s_all[segment_start_idx+col-1], nothing) for col in 1:stages]
            for stage = 2:stages
                fxSum = zeros(NonlinearExpr, nStates)
                for col = 1:stages
                    fxSum += tableau.a[stage, col] * F[col]
                end
                @constraint(model, X_raw[segment_start_idx+stage-1, :] .== X_raw[segment_start_idx, :] .+ h_all[segment] .* fxSum ./ x_scale)
                f(X[segment_start_idx+stage-1, :], U[segment_start_idx+stage-1, :], s_all[segment_start_idx+stage-1], model)
            end
            segment_start_idx = segment_start_idx + stages - 1
        end
        
        # ===================================== end of original =====================================
        #initialize all states and controls (set on raw O(1) variables)
        for idx = eachindex(s_all)
            init_vals = initialization.states(s_all[idx])
            init_u = initialization.controls(s_all[idx])
            set_start_value.(X_raw[idx, :], init_vals ./ x_scale)
            set_start_value.(U_raw[idx, :], init_u ./ u_scale)
        end

        return [model, X, U, s_all, segment_edges]
    end
 
    function create_interpolation(X_vals, U_vals, s_all, segment_edges)
        number_of_segments = length(segment_edges) - 1
        stride = stages - 1  # points per segment (adjacent segments share boundary)
 
        function find_segment(s)
            for i in 1:number_of_segments
                s <= segment_edges[i+1] && return i
            end
            return number_of_segments
        end
 

        function state_interp(s)
            seg   = find_segment(s)
            h     = segment_edges[seg+1] - segment_edges[seg]
            τ_eval = (s - segment_edges[seg]) / h       # map to [0,1]
            i0    = (seg - 1) * stride + 1
            i_end = i0 + stages - 1
            x_seg = X_vals[i0:i_end, :]
            return [_bary_eval(τ_ref, w_bary, x_seg[:, k], τ_eval) for k in 1:nStates]
        end
 
        function control_interp(s)
            seg    = find_segment(s)
            h      = segment_edges[seg+1] - segment_edges[seg]
            τ_eval = (s - segment_edges[seg]) / h
            i0     = (seg - 1) * stride + 1
            i_end  = i0 + stages - 1
            u_seg  = U_vals[i0:i_end, :]
            return [_bary_eval(τ_ref, w_bary, u_seg[:, k], τ_eval) for k in 1:nControls]
        end
 
        return Result_interpolation(state_interp, control_interp, s_all, collect(Float64, segment_edges))
    end
 
    LobattoIIIAMethod = Collocation(
        createDynamicConstraints,
        create_interpolation,
        nothing,
        nothing,
        f
    )
    return LobattoIIIAMethod
end




function smooth_by_OCP(track::Track, r::Float64, ds::Float64,closedTrack::Bool)
    println("Smoothing Track")
    x_smpl = track.x
    y_smpl = track.y

    dx = diff(x_smpl)
    dy = diff(y_smpl)
    ds_tmp = sqrt.(dx.^2 + dy.^2)

    s_tmp = [0; cumsum(ds_tmp)]

    ####
    N = Int(round((s_tmp[end] - s_tmp[1]) / ds))
    ds = (s_tmp[end] - s_tmp[1]) / N
    s_traj = LinRange(s_tmp[1], s_tmp[end], Int(N))
    x_smpl = interp1(s_tmp, x_smpl, s_traj,"PCHIP") 
    y_smpl = interp1(s_tmp, y_smpl, s_traj,"PCHIP") 
    ###

    dx = diff(x_smpl)
    dy = diff(y_smpl)

    th_init = unwrap(atan.(dy, dx))
    th_init = [th_init[1];th_init]
    C_init = diff(th_init) ./ ds
    C_init = [C_init; C_init[end]]

    ####
    # Initialize the optimization problem
    model = JuMP.Model(Ipopt.Optimizer)
    function f(z::AbstractVector, u::AbstractVector, s::Real, model)
        return [u[1]; z[1]; cos(z[2]); sin(z[2])]
    end

    Lobattostage = 2
    nStates = 4
    nControls = 1

    X_init = hcat(C_init, th_init, x_smpl, y_smpl)
    U_init = reshape(diff(C_init), :, 1)
    initialization = make_result_interpolation(X_init, U_init, collect(s_traj), collect(s_traj[1:end-1]))

    # Adaptive method defaults to bounds +/- scale, so choose generous scales from initialization.
    x_scale = 2.0 .* max.(vec(maximum(abs.(X_init), dims=1)), [1e-2, 1.0, 1.0, 1.0])
    u_scale = [2.0 * max(maximum(abs.(U_init)), 1e-3)]

    lobotom = createLobattoIIIA_Adaptive(f, Lobattostage, model, nControls, nStates, track;
        x_scale=x_scale, u_scale=u_scale)

    segment_edges = collect(s_traj)
    xd = lobotom.createConstraints(segment_edges, initialization)
    Xall = xd[2]
    Uall = xd[3]
    s_all = xd[4]
    segment_edges = xd[5]

    Z = Xall
    u = Uall[1:end-1, 1]
    s_nodes = s_all
    
    z_C = Z[:, 1]
    z_th= Z[:, 2]
    z_x = Z[:, 3]
    z_y = Z[:, 4]

    # Objective function (minimize the final time)
    x_dev = ((z_x[1:end-1] - x_smpl[1:end-1]).^2 + (z_x[2:end] - x_smpl[2:end]).^2) / 2
    y_dev = ((z_y[1:end-1] - y_smpl[1:end-1]).^2 + (z_y[2:end] - y_smpl[2:end]).^2) / 2

    @objective(model, Min, sum(ds .* (r .* u.^2 .+ x_dev .+ y_dev)))

    # Dynamic constraints
    if closedTrack
        @constraint(model, z_x[end] == z_x[1])
        @constraint(model, z_y[end] == z_y[1])
        @constraint(model, z_th[end] == z_th[1] + 2 * pi * round((th_init[end] - th_init[1]) / 2 / pi))
        @constraint(model, z_C[end] == z_C[1])
    else
        @constraint(model, z_th[end] == th_init[end])
        @constraint(model, z_x[end] == x_smpl[end])
        @constraint(model, z_y[end] == y_smpl[end])
        @constraint(model, z_C[end] == C_init[end])
    end

    # Solve NLP
    #set_silent(model)
    optimize!(model)
    # Extract solution values
#    @infiltrate

    x_traj = value.(z_x)
    y_traj = value.(z_y)
    C_traj = value.(z_C)
    th_traj = value.(z_th)

    itp = lobotom.createInterpolator(value(Xall), value(Uall), s_all, segment_edges)

    track.x = x_traj
    track.y = y_traj
    track.curvature = C_traj
    track.theta = th_traj
    track.sampleDistances = s_nodes
    track.fcurve = s -> begin
        vals = itp.states(s)
        return (vals[1], vals[2], vals[3], vals[4])
    end

    # Stretch track parameters
    s_new_endpoints = [s_nodes[1], s_nodes[end]]
    stretch(v) = length(v) == length(s_tmp) ?
        interp1(s_tmp, Float64.(v), s_nodes) :
        interp1(s_new_endpoints, [Float64(v[1]), Float64(v[1])], s_nodes)
    track.widthR      = stretch(track.widthR)
    track.widthL      = stretch(track.widthL)
    track.rho         = stretch(track.rho)
    track.μ           = stretch(track.μ)
    track.inclination = stretch(track.inclination)

    return (x_traj, y_traj, C_traj, th_traj)

end


function figureEight(vis::Bool = false, ds::Float64 = 0.5)
    w_l = 5.0  # Width of the track [m]
    w_r = 5.0  # Width of the track [m]
    rho = 1.225
    μ = 1.0

    w_x = 50.0  # width of the track along x axis
    w_y = 30.0  # width of the track along y axis
    t = collect(0:0.01:2π)
    X = w_x * cos.(t)
    Y = w_y * sin.(2 * t)

    track = Track(
        [0.0],      # curvature
        fill(rho, length(X)),
        fill(μ, length(X)),
        [1.0],      # sampleDistances
        trackMapping,
        X,
        Y,
        [0.0],      # theta
        fill(w_r, length(X)),
        fill(w_l, length(X)),
        [0.0],      # inclination
        [0.0],      # slope
        s -> (0.0), # fcurve
        [0.0]       # s
    )

    smooth_factor = 1e2
    smooth_by_OCP(track, smooth_factor, ds, true)

    if vis == 1
        plotTrack(track)
    end

    return track
end

function carODE_path(car::Car, track::Track, k::Union{Int64,Float64},
    u::AbstractVector, x::AbstractVector,
    model::Union{JuMP.Model,Nothing}=nothing)
    #car.mapping(car,u,x)
    car.controlMapping(u)
    car.stateMapping(x)
    dzds = time2path(car, track, k, model) #time2path(s,instantCarParams,track,car)
    return dzds
end;



function time2path(car::Car, track::Track, k::Union{Int64,Float64}, model::Union{Nothing,JuMP.Model})
    #track.mapping(track,instantTrack,s)
    dxdt = car.carFunction(track, model)
    v_x = car.carParameters.velocity.value[1]
    v_y = car.carParameters.velocity.value[2]
    psi = car.carParameters.psi.value
    n = car.carParameters.n.value

    th = track.fcurve(k)[2]
    C = track.fcurve(k)[1]

    epsilon = psi - th

    Sf = (1 - n * C) ./ (v_x .* cos(epsilon) - v_y .* sin(epsilon))
    dndt = v_x .* sin(epsilon) + v_y .* cos(epsilon)

    dzds = [
        Sf .* dxdt[1];    #dvx/ds
        Sf .* dxdt[2];    #dvy/ds
        Sf .* dxdt[3];    #dpsi/ds
        Sf .* dxdt[4];    #ddpsi/ds
        Sf .* dndt;       #dn/ds
        Sf               #dt/ds 
        # this should be capable of accepting more states
    ]

    #Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
    return dzds
end

function find_optimal_trajectory_adaptive(exp::Experiment, segments::Int64, pol_order::Int64, variant)
    track = exp.track
    car = exp.car
    model = exp.model
    performSensitivity = wants_sensitivity(exp.solver)

    #get size of control and state vector
    nControls = Int64(car.carParameters.nControls.value)
    nStates = Int64(car.carParameters.nStates.value)
    x = nothing
    u = nothing
    s_all = nothing

    function f(x, u, s, model)
        dxds = carODE_path(car, track, s, u, x, model)
        return dxds
    end

    #initialize solution = drive PD controller along centerline, one controller for all cars right now, but controller should be put inside car structure in future
    initialization = initializeSolution_interpolation(car, track, Int64(round(track.sampleDistances[end] * 2)))
    # Create sampling points for track
    segment_edges = LinRange(track.sampleDistances[1], track.sampleDistances[end], segments + 1)

    #From track that has data in vectors: vector of X, vector of Y vector of curvature ... create track that uses functions to describe its states
    # the interpoation function should be from the smoothing, but i am using basic interpolation of julia, work for future
    track_interp = interpolate_track(track)

    #this is used  to create constraints on state n, which gets its constraints from track width
    n_state_index = findfirst(e -> e.name == "n", car.carParameters.state_descriptor)

    clear = 0
    iterations = 0
    while clear == 0 && iterations <= 50
        clear = 1
        iterations += 1
        #@infiltrate

        # Scaling of controls and states, scaling is derived from bounds setup when creating a car
        # These functions 
        x_scale = get_scales(car.carParameters.state_descriptor)
        u_scale = get_scales(car.carParameters.control_descriptor)
        x_lb_static, x_ub_static = get_bounds(car.carParameters.state_descriptor)
        u_lb_static, u_ub_static = get_bounds(car.carParameters.control_descriptor)

        # state "n" needs special constraints
        margin = car.chassis.track.value/2
        function x_lb_fun(s)
            v = copy(x_lb_static)
            if n_state_index !== nothing
                v[n_state_index] = -track_interp.widthR(s) + margin
            end
            return v
        end
        function x_ub_fun(s)
            v = copy(x_ub_static)
            if n_state_index !== nothing
                v[n_state_index] = track_interp.widthL(s) - margin
            end
            return v
        end
        u_lb_fun(s) = u_lb_static
        u_ub_fun(s) = u_ub_static

        RKadaptive = createLobattoIIIA_Adaptive(f, pol_order, model, nControls, nStates, track;
            x_scale=x_scale, u_scale=u_scale,
            x_lb=x_lb_fun, x_ub=x_ub_fun, u_lb=u_lb_fun, u_ub=u_ub_fun)
        println("creating constraints")

        # Add parameters for sensitivity analysis, makes the problem little bigger, because parameters of car are now variables that ipopt sees
        tunables = carParameter[]
        if performSensitivity
            (params, tunables) = setParameters(car, model)
            exp.params = params
        end
        xd = RKadaptive.createConstraints(segment_edges, initialization)
        X = xd[2]
        U = xd[3]
        s_all = xd[4]
        segment_edges = xd[5]

        # Boundary conditions are supplied by the discipline.
        apply_boundary_conditions!(exp.discipline, model, X, track)

        # Global (lap-integral) constraints from the experiment.
        for gc in exp.global_constraints
            apply_global!(gc, model, X, U, car, s_all)
        end

        @constraint(model, diff(X[:, 6]) .>= 0) #time goes forward
        control_reg = 0#1e-4 * sum((U[i, j] / u_scale[j])^2 for i in axes(U, 1), j in axes(U, 2))
        @objective(model, Min, X[end, 6] + control_reg)

        #@constraint(model,-50 .<= diff(U[1:end,1])./diff(X[:,6]) .<= 50) #constraint on controls derivative
        #@constraint(model,-50 .<= diff(U[1:end,2])./diff(X[:,6]) .<= 50) #constraint on controls derivative
        #@constraint(model,-0.5 .<= diff(U[1:end-1,3]./diff(X[:,6])) .<= 0.5) #constraint on controls derivative

        println("Variables: $(num_variables(model)), Constraints: $(num_constraints(model; count_variable_in_set_constraints=true))")

        optimize!(model)
        println("Termination: ", termination_status(model), " | Objective: ", objective_value(model))

        x = value.(X)
        u = value.(U)
        exp.model = model
        # The resulting control is not just a vector for each state, but it is a function of path with fancy interpolation depending on discretization method
        exp.optiResult = RKadaptive.createInterpolator(x, u, s_all, segment_edges)

        # Here you can use the result of previous optimisiation as warm start, sometimes works very good sometimes very bad
        #initialization = make_result_interpolation(x, u, s_all)
        if performSensitivity
            println("resetting parameters for sensitivity analysis")
            exp.params = params
            resetParameters(tunables)
        end

        # Call to Mesh refinment algorigth, if clear == 1 the solution does not need to be refind and optimization ends
        #(segment_edges, clear, segment_errors) = refineMesh(exp, segment_edges, s_all, pol_order)
        clear = 1
        if clear == 0
            #create clear model for next iteration
            model = build_model(exp.solver)
        end
#        println("Sum of all errors: ", sum(segment_errors))
    end
    out = Result(x, u[1:end-1, :], s_all)
    return out, exp.optiResult
end

function interpolate_track(track::Track)
    s_nodes = collect(Float64.(track.sampleDistances))
    N = length(s_nodes)
    N >= 2 || error("interpolate_track: need at least two sample distances, got $(N)")
    #@infiltrate
    x_itp       = interpolate((s_nodes,), collect(Float64.(track.x)),           Gridded(Linear()))
    y_itp       = interpolate((s_nodes,), collect(Float64.(track.y)),           Gridded(Linear()))
    heading_itp = interpolate((s_nodes,), collect(Float64.(track.theta)),       Gridded(Linear()))
    curv_itp    = interpolate((s_nodes,), collect(Float64.(track.curvature)),   Gridded(Linear()))
    widthL_itp  = interpolate((s_nodes,), collect(Float64.(track.widthL)),      Gridded(Linear()))
    widthR_itp  = interpolate((s_nodes,), collect(Float64.(track.widthR)),      Gridded(Linear()))
    rho_itp     = interpolate((s_nodes,), collect(Float64.(track.rho)),         Gridded(Linear()))
    mu_itp      = interpolate((s_nodes,), collect(Float64.(track.μ)),           Gridded(Linear()))
    incl_itp    = interpolate((s_nodes,), collect(Float64.(track.inclination)), Gridded(Linear()))

    #x_itp       = Spline1D(s_nodes, collect(Float64.(track.x));           k=3, bc="nearest")
    #y_itp       = Spline1D(s_nodes, collect(Float64.(track.y));           k=3, bc="nearest")
    #heading_itp = Spline1D(s_nodes, collect(Float64.(track.theta));       k=3, bc="nearest")
    #curv_itp    = Spline1D(s_nodes, collect(Float64.(track.curvature));   k=3, bc="nearest")
    #widthL_itp  = Spline1D(s_nodes, collect(Float64.(track.widthL));      k=3, bc="nearest")
    #widthR_itp  = Spline1D(s_nodes, collect(Float64.(track.widthR));      k=3, bc="nearest")
    #rho_itp     = Spline1D(s_nodes, collect(Float64.(track.rho));         k=3, bc="nearest")
    #mu_itp      = Spline1D(s_nodes, collect(Float64.(track.μ));           k=3, bc="nearest")
    #incl_itp    = Spline1D(s_nodes, collect(Float64.(track.inclination)); k=3, bc="nearest")
    # Compatibility fcurve(s) -> (curvature, theta, x, y) matching make_fcurve order.
    fcurve_compat = track.fcurve

    return Track_interpolated(
        s_nodes,
        x_itp,
        y_itp,
        heading_itp,
        curv_itp,
        widthL_itp,
        widthR_itp,
        rho_itp,
        mu_itp,
        incl_itp,
        fcurve_compat,
    )
end

function initializeSolution_interpolation(car::Car, track::Track, segments::Int64; vref=5.0)
    println("started initialization")
    max_steer = car.wheelAssemblies[1].maxAngle.value
    # Constrants of controller
    Kv = 12 * car.carParameters.mass.value / 280.0
    Kp, Kd = 1.0   , 2.0
    nControls = Int64(car.carParameters.nControls.value)

    function ctrl(s, x)
        # PD controller + feedfw on steering, D is derived from lateral speed of car, cars heading and tracks heading
        # P controller on speed
        th = track.fcurve(s)[2]
        ε = atan(sin(x[3] - th), cos(x[3] - th))  # wrap to [-π, π]
        n_dot = x[1] * sin(ε) + x[2] * cos(ε)
        torque = (vref - x[1]) * Kv
        # feedforward curvature + feedback on lateral offset and heading
        C = track.fcurve(s)[1]
        δ_ff = atan(C * car.chassis.wheelbase.value)
        steering = clamp(δ_ff - (x[5] * Kp + n_dot * Kd), -max_steer, max_steer)
        return torque, steering
    end

    x0 = [vref, 0.0, track.theta[1], 0.0, 0.0, 0.0]
    s_span = (track.sampleDistances[1], track.sampleDistances[end])
    s_save = LinRange(s_span..., segments)

    # Live debug plot, because sometimes the initilization gets stuck(if the controller onstant are not good for vehicle)
    labels = ["vx", "vy", "ψ", "ψ̇", "n", "t", "torque", "steering"]
    #fig = Figure(size=(1200, 800))
    #debug_axes = [Axis(fig[div(i-1, 3)+1, mod(i-1, 3)+1], title=labels[i]) for i in eachindex(labels)]
    #debug_obs = [Observable(Point2f[]) for _ in eachindex(labels)]
    for i in eachindex(labels)
    #    lines!(debug_axes[i], debug_obs[i]; linewidth=2)
    end
    #display(GLMakie.Screen(), fig)

    cb = FunctionCallingCallback(; funcat=range(s_span..., length=segments)) do x, s, _
        pct = round(100 * (s - s_span[1]) / (s_span[2] - s_span[1]), digits=0)
        print("\r  Initialization: $(Int(pct))%")
        torque, steering = ctrl(s, x)
        for i in 1:6
            push!(debug_obs[i][], Point2f(s, x[i])); notify(debug_obs[i]); autolimits!(debug_axes[i])
        end
        for (j, v) in enumerate((torque, steering))
            push!(debug_obs[6+j][], Point2f(s, v)); notify(debug_obs[6+j]); autolimits!(debug_axes[6+j])
        end
    end

    # This is the actual ODEProblem which is solved car+track+controller
    # Note that problem is not solved in time but in PATH
    prob = ODEProblem((du, x, _, s) -> begin
        torque, steering = ctrl(s, x)
        du .= carODE_path(car, track, s, [torque, steering, zeros(nControls - 2)...], x, nothing)
    end, x0, s_span)

    # Plotting during initialization is super slow...
    #sol = OrdinaryDiffEq.solve(prob, Rodas4(autodiff=AutoFiniteDiff()), saveat=s_save, reltol=1e-4, abstol=1e-4, callback=cb)
    sol = OrdinaryDiffEq.solve(prob, Rodas4(autodiff=AutoFiniteDiff()), saveat=s_save, reltol=1e-4, abstol=1e-4)
    println()

    x = hcat(sol.u...)'
    s = sol.t
    u = zeros(segments, nControls)
    for i in eachindex(s)
        u[i, 1], u[i, 2] = ctrl(s[i], view(x, i, :))
    end

    # Finally interpolate the states and controls obtained during initialization for track
    initialization = make_result_interpolation(x, u, s)

    #Plot initialization
    #fig_path = Figure()
    #ax_path = Axis(fig_path[1, 1], aspect=DataAspect(), title="Initialization")
    #plotCarPath_interpolated(track, initialization, ax_path)
    #display(GLMakie.Screen(), fig_path)

    return initialization
end;

# =============================================START OF EXECUTION===========================================================================
track = figureEight(); 
car = createTwintrack();

ipopt_attrs = Dict{String,Any}(
    "linear_solver"          => "mumps",
    "print_timing_statistics"=> "yes",
    #"hessian_approximation" => "limited-memory",
    "alpha_for_y" => "safer-min-dual-infeas"
)


exp = Experiment(
    car = car,
    track = track,
    discipline = Open(v_start=5.0),            
    solver = IpoptBackend(
        performSensitivity = false,
        attributes = ipopt_attrs,
    ),
    global_constraints = GlobalConstraint[],    
    analysis = AnalysisConfig(
        plot_path       = true,
        plot_controls   = true,
        plot_states     = true,
        plot_jacobian   = true,
        plot_hessian    = true,
        animate         = false,
        animation_speedup = 1.0,
        time_simulation = false,
        sensitivity     = false,
    ),
)

segments = Int64(round(track.sampleDistances[end] *2)) #control size of the problem by seelcting segments per meter, now it is 2 segments per meter
pol_order = 2 # order of aproximating polynomial
run_experiment!(exp, segments, pol_order; variant="Lobatto"); # cretes problem and solves it
nothing