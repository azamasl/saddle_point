using LinearAlgebra, Convex, Random, CPUTime, Plots

include("global_solver.jl")
include("sp_function.jl")

n, m = 1, 1
TYPE = "2D"

#Initail point:
ran_seed, radius = 423327,10
Random.seed!(ran_seed);
obj2D = sad_point2D_objective()
@show x0 = radius * (rand(1) - rand(1))   #-3<x0<3
@show y0 = radius * (rand(1) - rand(1))   #-3<y0<3
# x0 = [50.0]
# y0 = [-10.0]
# x0=[0.0]
# y0=[0.0]


max_tr = 1000.0 #maximum allowed Δ
tr = 1  #initial Δ
max_it = 100#0.5*1e3
prt = 0# 0 don't print grad norm at every iterations; 1, do it.
F_tol = 1e-8
eta = 0.001

#settings for EGM:
reset_in_pt = 0
dummy = 0


@time @CPUtime x_sol, y_sol, iter, nfs, val, ngx, ngy, kk =
    tr_dogleg(x0, y0, obj2D, "dummy", max_it, prt, F_tol, tr, max_tr, eta)
println("tr_doglec results:")
println("number of iterations  = $kk")
println("L = $val")
println("|∇xL| = $ngx")
println("|∇yL| = $ngy")
#"Checking if the solution is an extermal point or a saddlepoint"
@show x_sol
@show y_sol
@show obj2D.∇xL(x_sol,y_sol)
@show obj2D.∇yL(x_sol,y_sol)
∇F = obj2D.∇F(x_sol,y_sol)
@show ∇F
∇L = ∇F
∇L[2,1] = -∇F[2,1]
∇L[2,2] = -∇F[2,2]
@show ∇L
eigVal, eigVec = eigen(∇L)
@show  eigVal



plot(
    range(1, iter, step = 1),
    nfs,
    yscale = :log10,
    label = "global_alg",
    ylabel = "||F||",
    title = string(TYPE, " m = $m, n=$n"),
)
