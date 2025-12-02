
using Optim
using Infiltrator
using LinearAlgebra
n = collect(1:10)


f(v) = sum(v)


## objective function,its gradient and hessian
function f_grad!(g, x)
    #fill!(g, 1.0)

    g[1:Int(length(g)/2)] .= 0
    g[Int(length(g)/2+1):length(g)] .= 1
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
        ## v̇²-v⁴c² <= 0 curvature lateral acceleration constraint
        #bottom_part[i,i] = 2*x[i]
        bottom_part[i,i+offset] = 2*x[offset + i] 
    end

    J = vcat(top_part,bottom_part)



    ## c₁'(x) = [1 / 3, -1, 0, 0]
    #J[1, 1] = 1 / 3
    #J[1, 2] = -1
    ## c₂'(x) = [0, 0, -1 / x[3]^2, -1 / x[4]^2]
    #J[2, 3] = -1 / x[3]^2
    #J[2, 4] = -1 / x[4]^2
    return J
end


function c_hess!(H, x, λ)
    offset = length(n)
    # c₁''(x) = 0
    # c₂''(x) = Diag([0, 0, 2 / x[3]^3, 2 / x[4]^2])
    #H[3, 3] += λ[2] * 2 / x[3]^3
    #H[4, 4] += λ[2] * 2 / x[4]^3
    for i = 1:length(n)
        H[offset + i,offset + i] += λ[offset + i] * 2
    end
    return H
end



curvature = zeros(5)
curvature = [curvature;1:5]
ay_max = 3
ax_max = 10

x0 =[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0]


lc1 = ones(length(n),1)*-ax_max
lc2 = -ay_max/curvature
lc = [lx1;lx2']


uc1 = ones(length(n),1)*ax_max
uc2 = ay_max/curvature
uc = [uc1;uc2']

lx = ones(2*length(n),1)*-999999
ux = ones(2*length(n),1)*999999


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
