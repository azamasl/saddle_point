using LinearAlgebra, Convex, Random,CPUTime, Plots
#function call_solver(prob_type, solver, ran_seed)
    #prob_type=0: bilinear
    #prob_type=1: quadratic
    # solver:
    # 1: alg 2 from Higher Order paper
    # 2: alg2 with first-order expansion in (14)
    # 5: alg2 with first-order expansion in (14), and D(z,zhat) instead of D(z,zt)
    #ran_seed: random number generator seed
    #n,m: size of the problem, x \in R^n, y \in R^m
    #T: number of iterations in HO paper
    #stepsize: stepsize in solver 3
    ran_seed=343
    prob_type=1
    solver= 1
    #n,m = 10,223
    #n,m = 37,45
    n,m = 3,5
    T = 49
    stepsize = 0.01
    max_it = 0.5*1e3
    prt = 0 # don't print grad norm at every iterations; 1, do it.
    ###################################

    include("sp_function.jl")
    include("solvers.jl")
    include("HO_paper.jl")
    Random.seed!(ran_seed);
    A = randn(m,n)
    B = zeros(n,n)
    C = zeros(m,m)
    if prob_type ==1#quadratic
        println(" \n The Problem is convex-concave")
        B = random_PSD(n)
        C = random_PSD(m)
    else
        println(" \n The Problem is bilinear")
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
    z = [x0;y0]
    if solver==1
        println("\n############ alg 2 ################")
        println("T = $T")
        @time @CPUtime  z = alg2(z,1e-6, T, obj,sp,0)
    elseif solver==2
        println("\n############ alg 2 with first-order eq. (14) ################")
        println("T = $T")
        @time @CPUtime  z = alg2(z,1e-6, T, obj,sp,1)
    elseif solver==5
        println("\n############ alg 2 with first-order eq. (14) and D(z,zhat) instead of D(z,zt) ################")
        @time @CPUtime  z = alg2(z,1e-6, T, obj,sp,2)
    else
        println("THIS THE END!")
    end
    x_sol = z[1:n]
    y_sol = z[n+1:n+m]
    val =obj.L(x_sol,y_sol)
    ngx = LinearAlgebra.norm(obj.∇xL(x_sol,y_sol))
    ngy = LinearAlgebra.norm(obj.∇yL(x_sol,y_sol))
    println("f value at the found solution = $val")
    println("|∇xL| = $ngx")
    println("|∇yL| = $ngy")
    println("Uncomment the next lines to print the solution")
    # display(x_sol)
    # display(y_sol)
