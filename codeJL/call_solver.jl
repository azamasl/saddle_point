#!/usr/bin/env julia

# Calling an indivitual solver
using LinearAlgebra, Convex, Random,CPUTime, Plots
include("problem_settings.jl")
solver = 7
c1 = 1e-4#Armijo parameter
do_ls = 1# 1: do line search, 2: use constant stepsize
stepsize = 0.01
max_it = 10#0.5*1e3
prt = 1 # 0 don't print grad norm at every iterations; 1, do it.
tol =1e-16
reset_in_pt = 0#"Set this to 1 to get diffrenet initial points, every time."



lb=[]
if  solver==3
    @time @CPUtime  x_sol, y_sol,it= Resolvent_alg(x0,y0,obj,stepsize,sp)
elseif solver==4
    @time @CPUtime  x_sol, y_sol,it = Newton_method(x0,y0,obj,sp)
elseif solver==6
    @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= secantUpdate_alg(x0,y0,obj,stepsize,sp,max_it,prt,reset_in_pt,do_ls)
    lb = ["Secant"]
elseif solver==7
    @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= secantShermanWoodbury_alg(x0,y0,obj,stepsize,sp,max_it,prt,reset_in_pt,do_ls)
    lb = ["Secant_INV"]
elseif solver==5
    @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= H_secant_alg(x0,y0,obj,stepsize,sp,max_it,prt,reset_in_pt,do_ls)
    lb = ["H_Secant"]
elseif  solver==8
    @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= Broyden(x0,y0,obj,stepsize,sp,max_it,prt,reset_in_pt,do_ls)
    lb = ["Broyden"]
elseif solver ==9
    @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= EGM(x0,y0,obj,stepsize,sp,max_it,prt,reset_in_pt,do_ls)
    lb = ["EGM"]
end

println("L = $val")
println("|∇xL| = $ngx")
println("|∇yL| = $ngy")
println("Uncomment the next lines to see the solution")
# display(x_sol)
# display(y_sol)
println("number of iterations  = $it")
#nf = Vector[normFs]
plot(range(1,it,step=1), normFs,yscale = :log10,label = lb[1], title = "log-scale plot of ||F|| at every iteration")
#typeof(normFs)

# ylabel = "Norm F",
# yscale = :log10,
# ytickfont = font(10, "Courier"),
# title = "log-scale plot of ||F|| at every iteration",
# label = "",#||F||
# color = :blue,  marker = (:dot, 1, 0.02))
