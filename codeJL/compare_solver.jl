#!/usr/bin/env julia

#Comparing different solvers
using LinearAlgebra, Convex, Random,CPUTime, Plots
include("problem_settings.jl")

"compare 1: Sec&Broy, 2:Sec&EGM"
compId = 1
c1 = 1e-4#Armijo parameter
do_ls = 1# 1: do line search, 2: use constant stepsize
stepsize = 0.01
max_it = 500#0.5*1e3
prt = 0 # 0 don't print grad norm at every iterations; 1, do it.
tol =1e-16
reset_in_pt = 0#"Set this to 1 to get diffrenet initial points, every time."
dis=1#00000

#"The plot package seems to be unstable and I could not get plotting 3 functions to work."
fun_num = 2
funs=Array{Function}(undef,fun_num) # compare 2 functions
lb=[]
nfs=Array{Vector}(undef,fun_num)
if compId ==1
    funs[1]=secantUpdate_alg
    funs[2]=Broyden
    lb = ["Secant" "Broyden"]
elseif compId == 2
    funs[1]=secantUpdate_alg
    funs[2]=EGM
    lb = ["Secant" "EGM"]
elseif compId ==3
    funs[1]=EGM
    funs[2]=Broyden
    lb = ["EGM" "Broyden"]
else
end

println("This is a dummy print and shouldn't be printed. For some reason my first print starting from here is not showing up. ")
let
    V = [] # Can contain anything
    #V = Array{Float64, 2}[] # Vector of Array{Float64,2}
    cc=1
    for f in funs
        println("Function $f is being called")
        @time @CPUtime  x_sol, y_sol,it, nfs[cc], val, ngx, ngy= f(x0,y0,obj,stepsize,sp,max_it,prt,reset_in_pt,do_ls)
        println("L = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        println("number of iterations  = $it")
        cc = cc+1
    end
end
# plot( 1:100, sqrt.(1:100), labels="Square Root")
# plot!( 1:10, sqrt.(1:10), labels="Se Root" )
l1 = size(nfs[1],1)
l2 = size(nfs[2],1)
plot(range(1,l1,step=1), nfs[1],yscale = :log10, label = lb[1])
plot!(range(1,l2,step=1), nfs[2],yscale = :log10,label = lb[2], title = "log-scale plot of ||F|| at every iteration")
#savefig("sample_plots/Name.png")
# plt = plot(V,
# ylabel = "Norm F",
# yscale = :log10,
# ytickfont = font(10, "Courier"),
# title = "log-scale plot of ||F|| at every iteration",
# label = lb,
# color = [:orange :blue],  marker = ([:dot :d], 1, 0.02))
#color = [:orange :blue :green],  marker = ([:dot :d :hex], 1, 0.02))
#gui(plt)
