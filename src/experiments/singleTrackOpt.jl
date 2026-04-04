using SLapSim
using GLMakie
import MathOptInterface as MOI
using UnoSolver
using UnicodePlots
using DiffOpt
#dark theme detector for linux KDE with kde-cli-tools installed
detect = Sys.islinux()
if detect
    x = "kreadconfig6"
    option1 = "--key"
    option2 = "LookAndFeelPackage"
    try
        xd = read(`$x $option1 $option2`, String) # remember the backticks ``
        is_dark = occursin("dark", lowercase(xd))
        println("Is dark theme: ", is_dark)
        if is_dark
            set_theme!(theme_dark())
        else
            set_theme!(Makie.current_default_theme())
        end
    catch
        println("Could not detect dark theme")
    end
end

update_theme!(
    fontsize = 11,
    fonts = (; regular = "Latin Modern Roman", bold = "Latin Modern Roman Bold"),
    palette = (color = Makie.to_colormap(:tab10),),
    colormap = :turbo,
    Axis = (
        titlesize = 13,
        xlabelsize = 12,
        ylabelsize = 12,
        xticklabelsize = 10,
        yticklabelsize = 10,
    ),
    Legend = (
        labelsize = 11,
        titlesize = 12,
    ),
    Colorbar = (
        ticklabelsize = 10,
        labelsize = 11,
    ),
)

GLMakie.closeall()

problem = Problem_config(nothing, nothing, nothing, nothing,nothing)



#car = createSimplestSingleTrack()
car = createTwintrack(true)
#car = createBus()
problem.car = car
path = "tracks/FSCZ.kml"
#track = figureEight(true, 0.1)
#track = singleTurn(50.0,5.0,true)
#track = doubleTurn(true,0.1)
track = kml2track(path, false, true)
#track = doubleTurn(false,0.1)
#track = skidpad(false)
problem.track = track

#UnoSolver.Optimizer

model = DiffOpt.nonlinear_diff_model(Ipopt.Optimizer)



JuMP.set_optimizer_attribute(model, "max_iter", 3000)
JuMP.set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
for (k,v) in [
    ("alpha_for_y", "safer-min-dual-infeas"),
    ("recalc_y", "yes"),
    ("recalc_y_feas_tol", 1e-4),
    ("adaptive_mu_globalization", "kkt-error"),
    ("quality_function_balancing_term", "cubic"),
    ("mu_min", 1e-8),
    ("nlp_scaling_constr_target_gradient", 1.0),
    ("nlp_scaling_obj_target_gradient", 1.0),
    ("nlp_scaling_min_value", 1e-6),
    ("jacobian_regularization_value", 1e-6),
    ("bound_relax_factor", 0.0),
    ("mumps_pivtol", 1e-4),
    ("min_refinement_steps", 2),
    ("max_refinement_steps", 20),
    #("acceptable_tol", 1e-5),
    #("acceptable_iter", 5),
    ("acceptable_dual_inf_tol", 1e-1)
]
    JuMP.set_optimizer_attribute(model, k, v)
end

problem.model = model
#model = JuMP.Model(() -> UnoSolver.Optimizer(preset="ipopt"))
#optiResult = findOptimalTrajectory(track,car,model,sampleDistances,initialization)
segments = Int64(round(track.sampleDistances[end]/2))
pol_order = 2
#optiResult, optiResult_interp = find_optimal_trajectory2(problem,segments,pol_order,"Radau")
optiResult, optiResult_interp = find_optimal_trajectory_adaptive(problem, segments, pol_order, "Lobatto")
problem.optiResult = optiResult_interp

if 1 == 1
    fig = Figure()
    ax = Axis(fig[1, 1], aspect=DataAspect())
    #plotCarPath(track,optiResult_interp,ax)
    plotCarPath_interpolated(track, optiResult_interp, ax)
    println(ax)
    screen = display(GLMakie.Screen(), fig)
    # simulate in time feed forward using optimal controls

    getError([1, 5], problem)
    #sol = timeSimulation(car, optiResult, track)
    sol = timeSimulation_interpolated(car, optiResult_interp, track)
    lines!(ax, getindex.(sol.u, 5), getindex.(sol.u, 6), label="Simulated in time")
    axislegend(ax, position=:rt)

    #@infiltrate
    fig = nothing
    ax = nothing
    SLapSim.plotCarStates_interp(optiResult_interp, 0.1)
    #SLapSim.plotCarStates2(optiResult)

    #include("../dataAnalysis/jacobian.jl")
    #include("../dataAnalysis/hessian_test2.jl")
    #println("jacobian")
    #println(UnicodePlots.spy(jacobian))
    #println("hessian")
    #println(UnicodePlots.spy(H_star))
    #fig_jac = Figure()
    #ax_jac = Axis(fig_jac[1, 1], title="Jacobian")
    #spy!(ax_jac, SparseArrays.sparse(rotr90(jacobian)))
    #display(GLMakie.Screen(), fig_jac)

  #  error_itp = getErrors(problem)
  #  plot(error_itp(optiResult_interp.path))
  #  plotErrorsOnTrack2D(problem; itp=error_itp)
#
  #  # Plot all controls + error in a single square window (stacked rows)
  #  fig = Figure(size=(900, 900))
  #  # 2x2 grid of plots, each plot gets its own colorbar to the right
  #  plots = [
  #      (s -> optiResult_interp.controls(s)[1], "moment_front"),
  #      (s -> optiResult_interp.controls(s)[1], "moment_rear"),
  #      (s -> optiResult_interp.controls(s)[2], "steering"),
  #      (s -> error_itp(s), "error"),
  #  ]

    # prepare error interpolant and append as a callable



   # # Layout: 2 rows × 4 columns (each plot + its colorbar occupies 2 columns)
   # for i in 1:length(plots)
   #     row = div(i - 1, 2) + 1
   #     col = mod(i - 1, 2) + 1
   #     ax_col = 2 * col - 1
   #     cb_col = 2 * col
   #     p = plots[i]
   #     f = p[1]
   #     name = p[2]
   #     ax = Axis(fig[row, ax_col], aspect=DataAspect())
   #     _, plt = plot_on_path(problem, f, name; axis=ax)
   #     Colorbar(fig[row, cb_col], plt; label=name, width=25)
   # end
   # display(GLMakie.Screen(), fig)
   # lines(error_itp(optiResult_interp.path), axis=(title="chyba v zavislosti na poloze", yscale=log10))
   # sampling_density = get_sampling_density(optiResult_interp.path)
   # plot_on_path(problem,sampling_density,"sampling density")

    
    snapshots = snapshot_car(car, optiResult_interp, track)
    #plot_parameters(snapshots, car,    ["drivetrain.motors[1].torque" , "drivetrain.motors[2].torque" ,"drivetrain.motors[3].torque","drivetrain.motors[4].torque"],"wheelAssemblies[1].steeringAngle")
    
    sensitivityAnalysis(problem)
    try
    plot_parameters(snapshots, car,    ["drivetrain.motors[1].torque" , "drivetrain.motors[2].torque" ],"wheelAssemblies[1].steeringAngle",["drivetrain.motors[3].torque" , "drivetrain.motors[4].torque" ])
    catch
        plot_parameters(snapshots, car,    ["drivetrain.motors[1].torque" , "drivetrain.motors[2].torque" ],"wheelAssemblies[1].steeringAngle")
    end
    animateCarDual(track, optiResult_interp, car; speedup=1, view_radius= car.chassis.wheelbase.value*3,cam_offset=3.0, savepath="results/animation.mp4")
end