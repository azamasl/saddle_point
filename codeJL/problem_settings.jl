include("sp_function.jl")
#include("solvers.jl")

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

#############Last settings#################
    #ran_seed=15 # The plots in the draft
    # ran_seed=13335
    # prob_type=0
    # m,n = 40,50
###############################
    ran_seed=12#2588907
    "0: random bilinear, 1: random quadratic CC, 2: random ill-cond. quadratic CC "
    prob_type=2
    n,m =  500,400  #200,300
    "Reciprocal condition number"
    rec_cond = 1e-3
    dis=1#00000
    ###################################
    Random.seed!(ran_seed);
    TYPE=""
    A = randn(m,n)
    B = zeros(n,n)
    C = zeros(m,m)
    if prob_type ==1        #quadratic
        println("The Problem is convex-concave")
        B = random_PSD(n)
        C = random_PSD(m)
        TYPE="quadratic convex-concave,"
    elseif prob_type ==0   #bilinear
        println("The Problem is bilinear")
        TYPE = "bilinear,"
    else   # ill-condined quadratic
        B = random_PSD_cond(n, rec_cond)
        C = random_PSD_cond(m, rec_cond)
        A = zeros(m,n)#random_rectangle_cond(n,m, rec_cond)
        TYPE="ill-cond. quadratic convex-concave,"
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
