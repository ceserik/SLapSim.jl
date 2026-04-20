# Adapter that turns a dimensional optimal-control problem into a unitless one.
struct ScaledProblem{F,XLB,XUB,ULB,UUB,I}
    f::F
    x_lb::XLB
    x_ub::XUB
    u_lb::ULB
    u_ub::UUB
    init::I
    x_scale::Vector{Float64}
    u_scale::Vector{Float64}
end

function scaledProblem(f, x_lb, x_ub, u_lb, u_ub, init, x_scale, u_scale)
    x_scale = collect(Float64, x_scale)
    u_scale = collect(Float64, u_scale)

    f_s = (y, v, s, model) -> f(y .* x_scale, v .* u_scale, s, model) ./ x_scale

    init_s = (
        states   = s -> init.states(s)   ./ x_scale,
        controls = s -> init.controls(s) ./ u_scale,
    )

    return ScaledProblem(
        f_s,
        s -> x_lb(s) ./ x_scale,
        s -> x_ub(s) ./ x_scale,
        s -> u_lb(s) ./ u_scale,
        s -> u_ub(s) ./ u_scale,
        init_s,
        x_scale,
        u_scale,
    )
end

unscale_x(sp::ScaledProblem, X_s) = X_s .* sp.x_scale'
unscale_u(sp::ScaledProblem, U_s) = U_s .* sp.u_scale'

# Print dimensional bounds + scales at a representative path position (default s=0).
function printBounds(sp::ScaledProblem, s::Real=0.0)
    xlb = sp.x_lb(s) .* sp.x_scale;  xub = sp.x_ub(s) .* sp.x_scale
    ulb = sp.u_lb(s) .* sp.u_scale;  uub = sp.u_ub(s) .* sp.u_scale
    r(v) = round(v; digits=3)
    println("bounds @ s=$(r(s))")
    for j in eachindex(xlb)
        println("  x[$j] [$(r(xlb[j])), $(r(xub[j]))]  scale=$(r(sp.x_scale[j]))")
    end
    for j in eachindex(ulb)
        println("  u[$j] [$(r(ulb[j])), $(r(uub[j]))]  scale=$(r(sp.u_scale[j]))")
    end
end

# Tighten the [-1,1] variable box on RK vars with path-dependent physical bounds.
function applyBounds!(sp::ScaledProblem, X_s, U_s, s_all)
    for i = eachindex(s_all)
        xlb = sp.x_lb(s_all[i]); xub = sp.x_ub(s_all[i])
        ulb = sp.u_lb(s_all[i]); uub = sp.u_ub(s_all[i])
        for j = axes(X_s, 2)
            set_lower_bound(X_s[i, j], xlb[j])
            set_upper_bound(X_s[i, j], xub[j])
        end
        for j = axes(U_s, 2)
            set_lower_bound(U_s[i, j], ulb[j])
            set_upper_bound(U_s[i, j], uub[j])
        end
    end
    return nothing
end
