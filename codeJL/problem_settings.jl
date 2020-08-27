include("sp_function.jl")
include("solvers.jl")

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
    #ran_seed: random number generator seed
    #n,m: size of the problem, x \in R^n, y \in R^m
    #stepsize: stepsize in solver 3
    #tol: gradient tolerance
    ran_seed=937#15
    prob_type=0
    m,n = 89,90
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
    x0 = dis*randn(n)
    y0 = dis*randn(m)