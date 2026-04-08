function initializeSolution_interpolation(car::Car, track::Track, segments::Int64; vref=5.0)
    println("started initialization")
    max_steer = car.wheelAssemblies[1].maxAngle.value
    Kv = 3 * car.carParameters.mass.value / 280.0
    Kp, Kd = 3.0   , 2.0
    nControls = Int64(car.carParameters.nControls.value)

    function ctrl(s, x)
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

    # Live debug plot
    labels = ["vx", "vy", "ψ", "ψ̇", "n", "t", "torque", "steering"]
    fig = Figure(size=(1200, 800))
    debug_axes = [Axis(fig[div(i-1, 3)+1, mod(i-1, 3)+1], title=labels[i]) for i in eachindex(labels)]
    debug_obs = [Observable(Point2f[]) for _ in eachindex(labels)]
    for i in eachindex(labels)
        lines!(debug_axes[i], debug_obs[i]; linewidth=2)
    end
    display(GLMakie.Screen(), fig)

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

    prob = ODEProblem((du, x, _, s) -> begin
        torque, steering = ctrl(s, x)
        du .= carODE_path(car, track, s, [torque, steering, zeros(nControls - 2)...], x, nothing)
    end, x0, s_span)

    sol = OrdinaryDiffEq.solve(prob, Rodas4(autodiff=AutoFiniteDiff()), saveat=s_save, reltol=1e-4, abstol=1e-4, callback=cb)
    #sol = OrdinaryDiffEq.solve(prob, Rodas4(autodiff=AutoFiniteDiff()), saveat=s_save, reltol=1e-4, abstol=1e-4)
    println()

    x = hcat(sol.u...)'
    s = sol.t
    u = zeros(segments, nControls)
    for i in eachindex(s)
        u[i, 1], u[i, 2] = ctrl(s[i], view(x, i, :))
    end

    initialization = make_result_interpolation(x, u, s)

    fig_path = Figure()
    ax_path = Axis(fig_path[1, 1], aspect=DataAspect(), title="Initialization")
    plotCarPath_interpolated(track, initialization, ax_path)
    display(GLMakie.Screen(), fig_path)

    return initialization
end;
