using Optim
using Infiltrator
using LinearAlgebra
using Plots
#using GLMakie
n = collect(1:10)

offset = length(n)  
f(v) = -sum(v[offset+1:end])


## objective function,its gradient and hessian
function f_grad!(g, x)

    g[1:offset] .= 0.0
    g[offset+1:length(g)] .= -1.0
    return g
end

function f_hess!(h, x)
    fill!(h, 0.0)
    
    return h
end



# constraint function, its gradient and hessian
function c(y, x)
    offset = length(n)
    y[1:offset] .= x[1:offset]
    for i = 1:offset
        y[offset + i] = x[offset + i]^2
    end

    for i = 1:offset-1
        y[2*offset + i] = x[offset + i]+ x[i] - x[offset + i + 1] 
    end

    return y
end

function c_jac!(J, x)
    offset = length(n)
    fill!(J, 0.0)
    for i = 1:offset
        J[i, i] = 1.0
    end
    for i = 1:offset
        J[offset + i, offset + i] = 2 * x[offset + i]
    end


    for i = 1:offset-1
        J[offset + offset + i,  i] = 1
        J[offset + offset + i,  offset +i] = 1
        J[offset + offset + i,  offset +i+1] = -1
    end
    
    return J

end

function c_hess!(H, x, λ)
    offset = length(n)
    
    for i = 1:offset
        H[offset + i, offset + i] = 2 * λ[offset + i]
    end
    return H
end



curvature = zeros(5).+0.0001
curvature = [curvature;1;0.001;0.001;2;0.001]

ay_max = 20
ax_max = 10

x0 =[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0]


lc1 = ones(length(n),1)*-ax_max
lc2 = -ay_max./curvature
lc = vec([lc1;lc2;zeros(length(n)-1)])


uc1 = ones(length(n),1)*ax_max
uc2 = ay_max./curvature
uc = vec([uc1;uc2;zeros(length(n)-1)])

lx = vec(ones(2*length(n),1)*-9999)
ux = vec(ones(2*length(n),1)*9999)


df = TwiceDifferentiable(f, f_grad!, f_hess!, x0)
dfc = TwiceDifferentiableConstraints(
    c, 
    c_jac!, 
    c_hess!, 
    lx, 
    ux, 
    lc, 
    uc,
)
res = optimize(df, dfc, x0, IPNewton())
xopt = Optim.minimizer(res)


# Plot solution: speeds vs max allowed, accelerations vs bounds
accel = xopt[1:offset]
speed = xopt[offset+1:end]
smax  = min.(100.0,sqrt.(ay_max ./ curvature))

plt = plot(layout=(2,1), size=(800,600))

plot!(plt[1], 1:offset, speed; label="speed", marker=:circle, xlabel="index", ylabel="speed")
plot!(plt[1], 1:offset, smax;  label="speed limit √(ay_max/curv)", ls=:dash)

plot!(plt[2], 1:offset, accel; label="accel", marker=:circle, xlabel="index", ylabel="accel")
hline!(plt[2], [ ax_max, -ax_max]; label=["+ax_max" "-ax_max"], ls=:dash)

display(plt)
savefig(plt, "solution_plots.png")
