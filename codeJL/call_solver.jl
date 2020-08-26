#!/usr/bin/env julia
using LinearAlgebra, Convex, Random,CPUTime, Plots

include("sp_function.jl")
include("solvers.jl")
include("HO_paper.jl")
#function call_solver(prob_type, solver, ran_seed)
    #prob_type=0: bilinear
    #prob_type=1: quadratic
    # solver:
    # 3: without extra-gradient step, i.e., z_{t+1} = z_t- s(s\nabla F(z_t)+I)^{-1} F(z_t), where s is the step-size
    # 4: Newton method: z_{t+1} = z_t- s(\nabla F(z_t))^{+} F(z_t), where + is the pseudo-inverse
    # 6: secant update alg
    # 7: secant update with extra gradient alg
    # 8: Broyden method
    # 9: EGM
    # 10: algs: it will run multiple solvers: 6,8(good),9
    #ran_seed: random number generator seed
    #n,m: size of the problem, x \in R^n, y \in R^m
    #stepsize: stepsize in solver 3
    #tol: gradient tolerance
    ran_seed=15
    "Set this to 1 to get diffrenet initial points, every time."
    reset_in_pt = 0

    prob_type=1
    n,m = 35,56

    M_opt = 1
    "If you set solver=10 it would compare two methods, See compId."
    solver= 10
    "1: Sec&Broy, 2:Sec&EGM"
    compId = 2

    stepsize = 0.01
    max_it = 500#0.5*1e3
    prt = 1 # 0 don't print grad norm at every iterations; 1, do it.
    tol =1e-16


    ###################################
    Random.seed!(ran_seed);

    A = randn(m,n)
    B = zeros(n,n)
    C = zeros(m,m)
    if prob_type ==1#quadratic
        println("The Problem is convex-concave")
        B = random_PSD(n)
        C = random_PSD(m)
    else
        println("The Problem is bilinear")
    end
    xstar = randn(n)
    ystar = randn(m)
    println("Uncomment the next lines to see the x^* and y^* ")
    #println("\n xstar and ystar:")
    #display(xstar)
    #display(ystar)
    sp = Saddle_Point(B,A,C,xstar,ystar)
    obj = saddle_point_objective(sp)
    x0 = randn(n)
    y0 = randn(m)

    if  solver==3
        @time @CPUtime  x_sol, y_sol,it= Resolvent_alg(x0,y0,obj,stepsize,sp)
    elseif solver==4
        @time @CPUtime  x_sol, y_sol,it = Newton_method(x0,y0,obj,sp)
    elseif solver==6
        @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= secantUpdate_alg(x0,y0,obj,stepsize,sp,M_opt,max_it,prt, reset_in_pt)
    elseif solver==7
        @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= secantUpdate_EGM(x0,y0,obj,stepsize,sp,max_it)
    elseif  solver==8
        @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= Broyden(x0,y0,obj,stepsize,sp,B_opt,max_it,prt,reset_in_pt)
    elseif solver ==9
        @time @CPUtime  x_sol, y_sol,it, normFs, val, ngx, ngy= EGM(x0,y0,obj,stepsize,max_it,prt,reset_in_pt)
    else
    end

    ###################print the final answer and plot
    if solver < 10
        println("L = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        println("Uncomment the next lines to see the solution")
        # display(x_sol)
        # display(y_sol)
        println("number of iterations  = $it")
        #nf = Vector[normFs]

        #typeof(normFs)

        # ylabel = "Norm F",
        # yscale = :log10,
        # ytickfont = font(10, "Courier"),
        # title = "log-scale plot of ||F|| at every iteration",
        # label = "",#||F||
        # color = :blue,  marker = (:dot, 1, 0.02))
end


###########################Comparing different algs
if solver==10
    "The plot package seems to be unstable and I could not get plotting 3 functions to work."
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
            @time @CPUtime  x_sol, y_sol,it, nfs[cc], val, ngx, ngy= f(x0,y0,obj,stepsize,sp,M_opt,max_it,prt,reset_in_pt)
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

end
#end
