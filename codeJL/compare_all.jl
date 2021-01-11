using LinearAlgebra, Convex, Random,CPUTime, JLD2
include("problem_settings.jl")
include("global_solver.jl")
include("solvers.jl")

"local alg settings"
c1 = 1e-4#Armijo parameter
do_ls = 1
#stepsize = 0.01
stepsize = 0.09#The fixed step_size, only used when do_ls = 0
"global alg settings"
tr = 1.0  #initial Δ
max_tr = 100.0 #maximum allowed Δ
eta = 0.001

"general settings"
fun_num = 5
max_it = 1000#.5*1e3
prt = 0 # 0 don't print grad norm at every iterations; 1, do it.
F_tol =1e-8
reset_in_pt = 0
dummy=0

nfs=Array{Vector}(undef,fun_num)
lb = ["Tr-Secant" "EGM" "Broyden" "Ls-Secant" "NoLs-Secant"]


@time @CPUtime  x_sol, y_sol,iter, nfs[1], val, ngx, ngy,kk = tr_dogleg(x0,y0, obj,sp, max_it, prt, F_tol, tr, max_tr, eta)
println("tr_doglec results:")
println("L = $val")
println("|∇xL| = $ngx")
println("|∇yL| = $ngy")
println("number of iterations  = $kk")

@time @CPUtime  x_sol, y_sol,iter, nfs[2], val, ngx, ngy = EGM(x0,y0,obj,dummy,sp,max_it,prt,reset_in_pt,dummy,F_tol)
println("EGM results:")
println("L = $val")
println("|∇xL| = $ngx")
println("|∇yL| = $ngy")
println("number of iterations  = $iter")

@time @CPUtime  x_sol, y_sol,iter, nfs[3], val, ngx, ngy = Broyden(x0,y0,obj,dummy,sp,max_it,prt,reset_in_pt,do_ls,F_tol)
println("Broyden results:")
println("L = $val")
println("|∇xL| = $ngx")
println("|∇yL| = $ngy")
println("number of iterations  = $iter")

@time @CPUtime  x_sol, y_sol,iter, nfs[4], val, ngx, ngy = secant_inv(x0,y0,obj,dummy,sp,max_it,prt,reset_in_pt,do_ls,F_tol)
println("local secant results:")
println("L = $val")
println("|∇xL| = $ngx")
println("|∇yL| = $ngy")
println("number of iterations  = $iter")

do_ls =0
@time @CPUtime  x_sol, y_sol,iter, nfs[5], val, ngx, ngy = secant_inv(x0,y0,obj,stepsize,sp,max_it,prt,reset_in_pt,do_ls,F_tol)
println("local secant results:")
println("L = $val")
println("|∇xL| = $ngx")
println("|∇yL| = $ngy")
println("number of iterations  = $iter")

#filname = "ws$prob_type.jld2"
@save "ws$prob_type.jld2" nfs lb TYPE prob_type ran_seed max_it F_tol stepsize
