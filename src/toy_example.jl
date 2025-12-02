
using Optim
using Infiltrator
using LinearAlgebra
n = collect(1:10)

offset = length(n)  
f(v) = sum(v[offset+1:end])


## objective function,its gradient and hessian
function f_grad!(g, x)
    #fill!(g, 1.0)

    g[1:offset] .= 0.0
    g[offset+1:length(g)] .= 1.0
    return g
end

function f_hess!(h, x)
    fill!(h, 0.0)
    
    return h
end



# constraint function, its gradient and hessian
function c(y, x)

    offset = length(n)
    # 10 <= v̇ <= 10
    y[1:offset] = x[1:offset]
    
    for i = 1:offset
        #@infiltrate
        ## -ay <= v̇² <= ay curvature lateral acceleration constraint
        y[offset+i] = x[offset + i]^2 
    end

    # y ma rozmer 2n
    # x ma rozmer 2n
    return y
end

 function c_jac!(J, x)
   # fill!(J, 0.0)
    offset = length(n)
    top_part = hcat(I(offset),zeros(offset,offset))

    bottom_part = zeros(10,20)
    for i = 1:offset
        #@infiltrate
        #bottom_part[i,i] = 2*x[i]
        bottom_part[i,i+offset] = 2*x[offset + i] 
    end

    J = vcat(top_part,bottom_part)

    return J
end


function c_hess!(H, x, λ)
    offset = length(n)

    for i = 1:length(n)
        H[offset + i,offset + i] += λ[offset + i] * 2
    end
    return H
end



curvature = zeros(5).+0.0001
curvature = [curvature;1:5]/5
ay_max = 20
ax_max = 10

x0 =[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0]


lc1 = ones(length(n),1)*-ax_max
lc2 = -ay_max./curvature
lc = vec([lc1;lc2])


uc1 = ones(length(n),1)*ax_max
uc2 = ay_max./curvature
uc = vec([uc1;uc2])

lx = vec(ones(2*length(n),1)*-999999)
ux = vec(ones(2*length(n),1)*999999)


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
optimize(df, dfc, x0, IPNewton())
