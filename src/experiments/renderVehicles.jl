include("transcription_benchmark.jl")
using GLMakie, CairoMakie

function _unique_car_factories(cases=_default_cases())
    seen = Set{String}()
    out = Tuple{String,Function}[]
    for case in cases
        for (name, fn) in case.cars
            if !(name in seen)
                push!(seen, name)
                push!(out, (name, fn))
            end
        end
    end
    return out
end

"""
    render_vehicle_matrix(; output_dir, track_for_car, cars)

Render a static top-down image for every unique vehicle in `_default_cases()`.
Each car factory is invoked once with a placeholder track (default: figureEight),
then `drawCar!` is called at origin with ψ=0. PDFs are saved as
`<output_dir>/<car_name>.pdf`.

Returns a list of named tuples `(car, path, success, [error])`.
"""
function render_vehicle_matrix(; output_dir::String="sync/vehicles",
    track_for_car=() -> figureEight(false, 0.1), cars=nothing)

    mkpath(output_dir)
    track = track_for_car()
    factories = _unique_car_factories()
    cars !== nothing && (factories = filter(c -> c[1] in cars, factories))

    results = []
    for (name, car_fn) in factories
        path = joinpath(output_dir, "$(name).pdf")
        entry = try
            car = car_fn(track)
            wb = car.chassis.wheelbase.value
            tw = car.chassis.track.value
            margin = max(wb, tw) * 1.2

            fig = Figure(size=(900, 700))
            ax = Axis(fig[1, 1], aspect=DataAspect(), title=name,
                xlabel="x [m]", ylabel="y [m]")

            chassis_obs = car.chassis.setupObservables(ax)
            wa_obs = [wa.setupObservables(ax, car.drivetrain.tires[i])
                      for (i, wa) in enumerate(car.wheelAssemblies)]
            aero_obs = car.aero.setupObservables(ax)

            ψ = π / 2  # car heading +y (up)
            Fz_static = car.chassis.mass.value * 9.81 / length(car.wheelAssemblies)
            car.chassis.updateObservables(chassis_obs, 0.0, 0.0, ψ)
            for (j, wa) in enumerate(car.wheelAssemblies)
                # Seed tire Fz so the friction ellipse has non-zero radius.
                car.drivetrain.tires[j].forces.value[3] = Fz_static
                wa.updateObservables(wa_obs[j], car.drivetrain.tires[j],
                    0.0, 0.0, ψ, Fz_static)
            end
            car.aero.updateObservables(aero_obs, 0.0, 0.0, ψ, wb, tw)

            xlims!(ax, -margin, margin)
            ylims!(ax, -margin, margin)
            CairoMakie.save(path, fig)
            (car=name, path=path, success=true)
        catch err
            (car=name, path=path, success=false, error=string(err))
        end
        push!(results, entry)
        @info "render done" car=entry.car success=entry.success path=entry.path
    end
    return results
end

render_vehicle_matrix()
