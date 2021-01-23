using CPUTime, JLD2, LinearAlgebra, Convex, Random
include("problem_settings.jl")
include("global_solver.jl")
include("solvers.jl")

#@btime

@time @CPUtime  x_sol[1], y_sol[1],iter[1], nfs[1], val[1], ng[1] = tr_dogleg(x0,y0, obj,sp, max_it, prt, F_tol, tr, max_tr, eta)

@time @CPUtime  x_sol[2], y_sol[2],iter[2], nfs[2], val[2], ng[2] = EGM(x0,y0,obj,dummy,sp,max_it,prt,dummy,F_tol)

@time @CPUtime  x_sol[3], y_sol[3],iter[3], nfs[3], val[3], ng[3] = Broyden(x0,y0,obj,dummy,sp,max_it,prt,do_ls,F_tol)

@time @CPUtime  x_sol[4], y_sol[4],iter[4], nfs[4], val[4], ng[4] = secant_inv(x0,y0,obj,dummy,sp,max_it,prt,do_ls,F_tol)

@time @CPUtime  x_sol[5], y_sol[5],iter[5], nfs[5], val[5], ng[5] = secant_inv(x0,y0,obj,stepsize,sp,max_it,prt,do_ls=0,F_tol)

#x_sol[1], y_sol[1],iter[1], nfs[1], val[1], ng[1] = tr_dogleg(x0,y0, obj,sp, max_it, prt, F_tol, tr, max_tr, eta)
#x_sol[4], y_sol[4],iter[4], nfs[4], val[4], ng[4] = secant_inv(x0,y0,obj,dummy,sp,max_it,prt,do_ls,F_tol)

#filname = "ws$prob_type.jld2"
@save "ws$prob_type.jld2" nfs lb TYPE prob_type ran_seed max_it F_tol stepsize x_sol y_sol iter val ng
