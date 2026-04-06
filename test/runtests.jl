using SLapSim
using Test

@testset "SLapSim.jl" begin
    @testset "Track_interpolated converter" begin
        # Build a minimal Track by hand with varying widths and known fields.
        # Avoid running the OCP smoothing path so the test stays fast/deterministic.
        N = 11
        s_nodes = collect(range(0.0, 10.0; length = N))
        x_vals  = collect(range(0.0, 20.0; length = N))
        y_vals  = collect(range(5.0,  6.0; length = N))
        th_vals = collect(range(0.0,  1.0; length = N))
        C_vals  = collect(range(-0.1, 0.2; length = N))
        wL_vals = collect(range(1.0,  3.0; length = N))   # varying left width
        wR_vals = collect(range(2.5,  1.5; length = N))   # varying right width
        rho_vals = fill(1.225, N)
        mu_vals  = fill(1.0,   N)
        incl_vals = zeros(N)

        track = Track(
            C_vals,                       # curvature
            rho_vals,                     # rho
            mu_vals,                      # μ
            s_nodes,                      # sampleDistances
            (a, b, s) -> nothing,         # mapping (unused here)
            x_vals, y_vals, th_vals,
            wR_vals, wL_vals,
            incl_vals,
            zeros(N),                     # slope
            s -> (0.0, 0.0, 0.0, 0.0),    # fcurve placeholder
            s_nodes                       # s
        )

        ti = interpolate_track(track)

        # converter correctness at every node: each callable must reproduce the
        # underlying value at the sample point exactly.
        @test ti.s_nodes == s_nodes
        for (i, s) in enumerate(s_nodes)
            @test ti.x(s)         ≈ x_vals[i]
            @test ti.y(s)         ≈ y_vals[i]
            @test ti.heading(s)   ≈ th_vals[i]
            @test ti.curvature(s) ≈ C_vals[i]
            @test ti.widthL(s)    ≈ wL_vals[i]
            @test ti.widthR(s)    ≈ wR_vals[i]
            @test ti.air_density(s) ≈ rho_vals[i]
            @test ti.friction_c(s)  ≈ mu_vals[i]
        end

        # Compatibility fcurve(s) -> (curvature, theta, x, y) at a node.
        s_mid = s_nodes[5]
        c, th, xm, ym = ti.fcurve(s_mid)
        @test c  ≈ C_vals[5]
        @test th ≈ th_vals[5]
        @test xm ≈ x_vals[5]
        @test ym ≈ y_vals[5]

        # Linear interpolation between nodes (basic sanity).
        s_between = (s_nodes[3] + s_nodes[4]) / 2
        @test ti.widthL(s_between) ≈ (wL_vals[3] + wL_vals[4]) / 2
        @test ti.widthR(s_between) ≈ (wR_vals[3] + wR_vals[4]) / 2
    end

    @testset "Envelope-based static n bounds" begin
        # Build a Track with varying widths; verify createTwintrack/formulaE2026
        # take the envelope (max) over the whole track instead of widthL[1].
        N = 7
        s_nodes = collect(range(0.0, 6.0; length = N))
        wL = [1.5, 2.0, 3.0, 4.0, 3.5, 2.5, 1.5]
        wR = [1.5, 1.7, 2.2, 2.8, 2.0, 1.6, 1.5]
        track = Track(
            zeros(N), fill(1.225, N), fill(1.0, N), s_nodes,
            (a, b, s) -> nothing,
            collect(range(0.0, 12.0; length = N)),
            zeros(N), zeros(N),
            wR, wL, zeros(N), zeros(N),
            s -> (0.0, 0.0, 0.0, 0.0), s_nodes
        )

        margin = 0.6
        for car in (createTwintrack(true, track), formulaE2026(track))
            n_entry = first(filter(e -> e.name == "n", car.carParameters.state_descriptor))
            @test n_entry.lb ≈ -maximum(wR) + margin
            @test n_entry.ub ≈  maximum(wL) - margin
        end
    end

    @testset "Track_interpolated bounds vary along s" begin
        # When width varies, the bounds derived from Track_interpolated must
        # change with s — this is what the adaptive solver consumes per node.
        N = 5
        s_nodes = collect(range(0.0, 4.0; length = N))
        wL = [2.0, 2.5, 3.5, 2.5, 2.0]
        wR = [2.0, 2.0, 2.0, 2.0, 2.0]
        track = Track(
            zeros(N), fill(1.225, N), fill(1.0, N), s_nodes,
            (a, b, s) -> nothing,
            zeros(N), zeros(N), zeros(N),
            wR, wL, zeros(N), zeros(N),
            s -> (0.0, 0.0, 0.0, 0.0), s_nodes
        )
        ti = interpolate_track(track)
        margin = 0.6

        n_lb = s -> -ti.widthR(s) + margin
        n_ub = s ->  ti.widthL(s) - margin

        # Tightest at the ends, loosest in the middle.
        @test n_ub(s_nodes[1]) ≈ wL[1] - margin
        @test n_ub(s_nodes[3]) ≈ wL[3] - margin
        @test n_ub(s_nodes[3]) > n_ub(s_nodes[1])
        @test n_lb(s_nodes[3]) ≈ -wR[3] + margin
    end
end
