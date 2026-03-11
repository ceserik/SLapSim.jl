using GLMakie

# Define the function and its gradient
f(x1, x2) = x1^2 + x2^2
∇f(x1, x2) = (2x1, 2x2)

# Constraint: g(x1, x2) = x1 + x2 + 3 = 0  →  x2 = -x1 - 3
g(x1, x2)  = x1 + x2 + 3
∇g(x1, x2) = (1.0, 1.0)   # constant gradient

# ── Analytical optimum (KKT: ∇f = λ∇g, g = 0) ───────────────────────────
# x1* = x2* = -3/2,  λ* = -3
x_opt, y_opt = -1.5, -1.5
λ_opt = -3.0

# Grid for contour/heatmap
n = 200
xs = LinRange(-5, 3, n)
ys = LinRange(-5, 3, n)
Z = [f(x, y) for y in ys, x in xs]

# Coarser grid for objective gradient arrows
ng = 12
xg = LinRange(-4.5, 2.5, ng)
yg = LinRange(-4.5, 2.5, ng)
X_pts = [x for y in yg, x in xg] |> vec
Y_pts = [y for y in yg, x in xg] |> vec
dX = [∇f(x, y)[1] for (x, y) in zip(X_pts, Y_pts)]
dY = [∇f(x, y)[2] for (x, y) in zip(X_pts, Y_pts)]

# Normalize objective gradient arrows for uniform length display
norms = sqrt.(dX .^ 2 .+ dY .^ 2)
scale = 0.25
dX_n = scale .* dX ./ norms
dY_n = scale .* dY ./ norms

# Filter out small red arrows near the optimum so they don't overlap the big arrows
clear_radius = 1.2
keep = sqrt.((X_pts .- x_opt).^2 .+ (Y_pts .- y_opt).^2) .> clear_radius
X_pts, Y_pts = X_pts[keep], Y_pts[keep]
dX_n,  dY_n  = dX_n[keep],  dY_n[keep]

# Constraint line: extend well beyond the visible axis limits
nc = 14
xc = LinRange(-7.0, 4.0, nc)
yc = @. -xc - 3.0
# Constraint gradient arrows (normalised, perpendicular to constraint line)
gx, gy = ∇g(0.0, 0.0)   # constant
gnorm  = sqrt(gx^2 + gy^2)
cgx = scale * gx / gnorm
cgy = scale * gy / gnorm

# Plot
fig = Figure(size = (750, 650))
ax = Axis(fig[1, 1];
    xlabel = "x₁", ylabel = "x₂",
    #title  = "f = x₁²+x₂², constraint x₁+x₂+3=0",
    aspect = DataAspect(),
    limits = (-4.0, 1.5, -4.0, 1.5))

# Filled contour (top view)
cf = contourf!(ax, xs, ys, Z; levels = 20, colormap = :viridis)
contour!(ax, xs, ys, Z; levels = 20, color = :white, linewidth = 0.4, alpha = 0.5)

# Objective gradient arrows (red)
arrows!(ax, X_pts, Y_pts, dX_n, dY_n;
    color     = :red,
    arrowsize = 8,
    linewidth = 1.5)

# Constraint line
lines!(ax, xc, yc; color = :orange, linewidth = 3)

# Constraint gradient arrows (orange, along the line)
arrows!(ax, collect(xc), collect(yc), fill(cgx, nc), fill(cgy, nc);
    color     = :orange,
    arrowsize = 10,
    linewidth = 2.0)

# ── Big arrows at the optimum ────────────────────────────────────────────
# Proportional lengths: |∇f| = |λ*| · |∇g| = 3 × |∇g|, opposite direction
unit = 0.55          # display length for one |∇g| unit (data coords)
d    = 1.0 / sqrt(2.0)   # equal component along each axis for (±1,±1) direction

# ∇g arrow (darkorange): direction (+1,+1)/√2, length = unit
arrows!(ax, [x_opt], [y_opt], [unit * d], [unit * d];
    color = :darkorange, arrowsize = 22, linewidth = 7.0)

# ∇f arrow (crimson): direction (-1,-1)/√2, length = 3*unit  (since λ* = -3)
arrows!(ax, [x_opt], [y_opt], [-3unit * d], [-3unit * d];
    color = :crimson, arrowsize = 22, linewidth = 7.0)

# Tick marks on ∇f at 1× and 2× unit — dividing it into 3 equal |∇g|-length pieces
perp = 0.13   # half-length of tick, perpendicular direction is (d, -d)
for k in 1:2
    tx = x_opt - k * unit * d
    ty = y_opt - k * unit * d
    linesegments!(ax,
        [tx - perp * d, tx + perp * d],
        [ty + perp * d, ty - perp * d];
        color = :white, linewidth = 2.5)
end

# Optimum point marker (drawn last so it sits on top)
scatter!(ax, [x_opt], [y_opt];
    color       = :white,
    strokecolor = :black,
    strokewidth = 2,
    markersize  = 14)

# Arrow labels
text!(ax, x_opt + unit * d + 0.05, y_opt + unit * d;
    text = "∇g", color = :darkorange, fontsize = 28, font = :bold)
text!(ax, x_opt - 3unit * d - 0.05, y_opt - 3.5unit * d - 0.22;
    text = "∇f", color = :crimson, fontsize = 28, font = :bold)

# λ annotation alongside the ∇f arrow (offset to the right of it)
text!(ax, x_opt - 1.5unit * d + 0.18, y_opt - 3unit * d + 0.08;
    text = "λ* = −3", color = :white, fontsize = 28, font = :bold)
text!(ax, x_opt - 1.5unit * d + 0.18, y_opt - 3unit * d - 0.18;
    text = "∇f = λ*∇g", color = :white, fontsize = 28)
# Legend
elem_f     = LineElement(color = :red,        linewidth = 2)
elem_g     = LineElement(color = :orange,     linewidth = 2)
elem_∇f    = LineElement(color = :crimson,    linewidth = 3)
elem_∇g    = LineElement(color = :darkorange, linewidth = 3)
#Legend(fig[1, 3],
#    [elem_f, elem_g, elem_∇f, elem_∇g],
#    ["∇f field", "constraint & ∇g field", "∇f at x*", "∇g at x*"])
#Colorbar(fig[1, 2], cf; label = "f(x₁, x₂)")

display(fig)
