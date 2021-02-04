using CPUTime, JLD2, LinearAlgebra, Convex, Random
include("problem_settings.jl")
include("global_solver.jl")
include("solvers.jl")
#include("global_solver_test.jl")

# x_sol, y_sol,iter, nfs, val, ng = tr_dogleg(x0,y0, obj,sp, max_it, prt, F_tol, tr, max_tr, eta)
# ret=fun_return(x_sol, y_sol, iter, nfs, val, ng)
# @save "output/ws_tr_$prob_type-$ts.jld2"  TYPE prob_type ran_seed max_it F_tol stepsize ret
#
# x_sol, y_sol,iter, nfs, val, ng = EGM(x0,y0,obj,dummy,sp,max_it,prt,dummy,F_tol)
# ret=fun_return(x_sol, y_sol, iter, nfs, val, ng)
# @save "output/ws_EGM_$prob_type-$ts.jld2"  TYPE prob_type ran_seed max_it F_tol stepsize ret
#
#
# x_sol, y_sol,iter, nfs, val, ng  = secant_inv(x0,y0,obj,dummy,sp,max_it,prt,do_ls,F_tol)
# ret=fun_return(x_sol, y_sol, iter, nfs, val, ng)
# @save "output/ws_lsSec_$prob_type-$ts.jld2" TYPE prob_type ran_seed max_it F_tol stepsize ret
#
#
# # x_sol, y_sol,iter, nfs, val, ng = tr_NW41(x0,y0, obj,sp, max_it, prt, F_tol, tr, max_tr, eta)
# # ret=fun_return(x_sol, y_sol, iter, nfs, val, ng)
# # @save "output/ws_tr_NW41_$prob_type-$ts.jld2" TYPE prob_type ran_seed max_it F_tol stepsize ret
#
#
x_sol, y_sol,iter, nfs, val, ng = Broyden(x0,y0,obj,dummy,sp,max_it,prt,do_ls,F_tol)
ret=fun_return(x_sol, y_sol, iter, nfs, val, ng)
@save "output/ws_BROYDEN_$prob_type-$ts.jld2"  TYPE prob_type ran_seed max_it F_tol stepsize ret
#
#
# x_sol, y_sol,iter, nfs, val, ng = secant_inv(x0,y0,obj,stepsize,sp,max_it, prt,      0,F_tol)
# ret=fun_return(x_sol, y_sol, iter, nfs, val, ng)
# @save "output/ws_noSec_$prob_type-$ts.jld2"  TYPE prob_type ran_seed max_it F_tol stepsize ret


x_sol, y_sol,iter, nfs, val, ng = Broyden(x0,y0,obj,stepsize,sp,max_it,prt,     0,F_tol)
ret=fun_return(x_sol, y_sol, iter, nfs, val, ng)
@save "output/ws_noBROYDEN_$prob_type-$ts.jld2"  TYPE prob_type ran_seed max_it F_tol stepsize ret
