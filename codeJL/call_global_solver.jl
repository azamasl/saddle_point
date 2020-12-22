using LinearAlgebra, Convex, Random,CPUTime, Plots
include("problem_settings.jl")
include("global_solver.jl")
max_tr = 100.0 #maximum allowed Δ
tr = 1.0  #initial Δ
max_it = 20#0.5*1e3
prt = 0 # 0 don't print grad norm at every iterations; 1, do it.
F_tol =1e-8
eta = 0.001

# lb=[]
# nfs=Array{Vector}(undef,fun_num)

x_sol, y_sol,iter, nfs, val, ngx, ngy,kk= tr_dogleg(x0,y0, obj,sp, max_it, prt, F_tol, tr, max_tr, eta)
println("############ Tr-Dogleg Method ###############")
println("L = $val")
println("|∇xL| = $ngx")
println("|∇yL| = $ngy")
println("Uncomment the next lines to see the solution")
# display(x_sol)
# display(y_sol)
println("number of total iterations  = $kk")


plot(range(1,iter,step=1), nfs, yscale = :log10,label = "global_alg", ylabel = "||F||", title = string(TYPE," m = $m, n=$n"))
#plot!(range(1,l2,step=1), nfs[2],yscale = :log10,label = lb[2], ylabel = "log ||F||", title = string(TYPE," m = $m, n=$n"))
#savefig("sample_plots/Name.png")
