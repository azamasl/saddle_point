#function call_solver(prob_type, solver, ran_seed)
    #prob_type=0: bilinear
    #prob_type=1: quadratic
    # solver:
    # 1: alg 2 from Higher Order paper
    # 2: alg2 with first-order expansion in (14)
    # 3: without extra-gradient step, i.e., z_{t+1} = z_t- s(s\nabla F(z_t)+I)^{-1} F(z_t), where s is the step-size
    # 4: Newton method: z_{t+1} = z_t- s(\nabla F(z_t))^{+} F(z_t), where + is the pseudo-inverse
    # 5: alg2 with first-order expansion in (14), and D(z,zhat) instead of D(z,zt)
    # 6: secant update alg
    # 7: secant update with extra gradient alg
    #ran_seed: random number generator seed
    #n,m: size of the problem, x \in R^n, y \in R^m
    #T: number of iterations in HO paper
    #stepsize: stepsize in solver 3
    prob_type=1
    solver=7
    ran_seed=3653
    #n,m = 10,223
    n,m = 102,123
    T = 49
    stepsize = 0.1
    ###################################
    using LinearAlgebra, Convex, Random,CPUTime
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

    if solver==4
        println("\n############ Newton Method ###############")
        @time @CPUtime x_sol, y_sol,it=Newton_method(x0,y0,obj,sp)
        val =obj.L(x_sol,y_sol)
        ngx = LinearAlgebra.norm(obj.∇xL(x_sol,y_sol))
        ngy = LinearAlgebra.norm(obj.∇yL(x_sol,y_sol))
        println("L = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        println("Uncomment the next lines to see the solution")
        # display(x_sol)
        # display(y_sol)
        println("number of iterations  = $it")

    elseif solver==1
        println("\n############ alg 2 ################")
        println("T = $T")
        z = [x0;y0]
        @time @CPUtime  z = alg2(z,1e-6, T, obj,sp,0)
        x_sol = z[1:n]
        y_sol = z[n+1:n+m]
        val =obj.L(x_sol,y_sol)
        ngx = LinearAlgebra.norm(obj.∇xL(x_sol,y_sol))
        ngy = LinearAlgebra.norm(obj.∇yL(x_sol,y_sol))
        println("f value at the found solution = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        println("Uncomment the next lines to see the solution")
        # display(x_sol)
        # display(y_sol)
    elseif solver==2
        println("\n############ alg 2 with first-order eq. (14) ################")
        z = [x0;y0]
        println("T = $T")
        @time @CPUtime  z = alg2(z,1e-6, T, obj,sp,1)
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
    elseif solver==5
        println("\n############ alg 2 with first-order eq. (14) and D(z,zhat) instead of D(z,zt) ################")
        z = [x0;y0]
        @time @CPUtime  z = alg2(z,1e-6, T, obj,sp,2)
        x_sol = z[1:n]
        y_sol = z[n+1:n+m]
        val =obj.L(x_sol,y_sol)
        ngx = LinearAlgebra.norm(obj.∇xL(x_sol,y_sol))
        ngy = LinearAlgebra.norm(obj.∇yL(x_sol,y_sol))
        println("f value at the found solution = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        println("Uncomment the next lines to print the solution")
    elseif  solver==6
        println("\n ############ Using Secant Update ############")
        @time @CPUtime  x_sol, y_sol,it= secantUpdate_alg(x0,y0,obj,0.01,sp)
        val =obj.L(x_sol,y_sol)
        ngx = LinearAlgebra.norm(obj.∇xL(x_sol,y_sol))
        ngy = LinearAlgebra.norm(obj.∇yL(x_sol,y_sol))
        println("L = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        println("Uncomment the next lines to see the solution")
        # display(x_sol)
        # display(y_sol)
        println("number of iterations  = $it")
    elseif  solver==7
        println("\n ############ Using Secant Update With Extra Gradient ############")
        @time @CPUtime  x_sol, y_sol,it= secantUpdate_extra_grad_alg(x0,y0,obj,0.01,sp)
        val =obj.L(x_sol,y_sol)
        ngx = LinearAlgebra.norm(obj.∇xL(x_sol,y_sol))
        ngy = LinearAlgebra.norm(obj.∇yL(x_sol,y_sol))
        println("L = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        println("Uncomment the next lines to see the solution")
        # display(x_sol)
        # display(y_sol)
        println("number of iterations  = $it")
    else # solver==3
        println("\n ############ Using Proximal Point ############")
        @time @CPUtime  x_sol, y_sol,it= Resolvent_alg(x0,y0,obj,stepsize,sp)
        val =obj.L(x_sol,y_sol)
        ngx = LinearAlgebra.norm(obj.∇xL(x_sol,y_sol))
        ngy = LinearAlgebra.norm(obj.∇yL(x_sol,y_sol))
        println("L = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        println("Uncomment the next lines to see the solution")
        # display(x_sol)
        # display(y_sol)
        println("number of iterations  = $it")
     end
#end
