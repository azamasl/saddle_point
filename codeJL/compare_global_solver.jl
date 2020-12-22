using LinearAlgebra, Convex, Random,CPUTime, Plots
include("problem_settings.jl")
include("global_solver.jl")
include("solvers.jl")
max_tr = 100.0 #maximum allowed Δ
tr = 1.0  #initial Δ
max_it = 700#0.5*1e3
prt = 0 # 0 don't print grad norm at every iterations; 1, do it.
F_tol =1e-8
eta = 0.001

#settings for EGM:
reset_in_pt = 0
dummy=0

fun_num = 2
funs=Array{Function}(undef,fun_num) # compare 2 functions
lb=[]
nfs=Array{Vector}(undef,fun_num)
funs[1]=tr_dogleg
funs[2]=EGM
lb = ["Tr-Secant" "EGM"]


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


l1 = size(nfs[1],1)
l2 = size(nfs[2],1)
plot(range(1,l1,step=1), nfs[1], yscale = :log10, label = lb[1])
plot!(range(1,l2,step=1), nfs[2], yscale = :log10,label = lb[2],xlabel = "Iteration", ylabel = "||F||", title = string(TYPE," m = $m, n=$n"))
#plot!(range(1,l2,step=1), nfs[2],yscale = :log10,label = lb[2], ylabel = "log ||F||", title = string(TYPE," m = $m, n=$n"))
savefig("sample_plots/globalName.png")
