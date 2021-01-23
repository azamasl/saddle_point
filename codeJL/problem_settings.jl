include("sp_function.jl")
#include("solvers.jl")

"general settings"
    ran_seed=102333#2588907
    "0: random bilinear, 1: random quadratic CC, 2: random ill-cond. quadratic CC "
    prob_type=0
    n,m =  120,100# 120,100#
    ts = m+n #total size
    "Reciprocal condition number"
    rec_cond = 1e-3
    dis=1
    fun_num = 5
    max_it = 5000#.5*1e3
    prt = 0 # 0 don't print grad norm at every iterations; 1, do it.
    F_tol =1e-8
    reset_in_pt = 0
    dummy=0

"local alg settings"
    c1 = 1e-4#Armijo parameter
    do_ls = 1
    #stepsize = 0.01
    stepsize = 0.09#The fixed step_size, only used when do_ls = 0

"global alg settings"
    tr = 0.1  #initial Δ
    max_tr = 0.5 #maximum allowed Δ
    eta = 0.01

lb = ["Tr-Secant" "EGM" "Broyden" "Ls-Secant" "NoLs-Secant"]

nfs=Array{Vector}(undef,fun_num)
x_sol=Array{Vector}(undef,fun_num)
y_sol=Array{Vector}(undef,fun_num)
iter=[]
val =[]
ng =[]

"Creating the random problme instance"
    Random.seed!(ran_seed);
    xstar = randn(n)
    ystar = randn(m)
    println("Uncomment the next lines to see the x^* and y^* ")
    #println("\n xstar and ystar:")
    #display(xstar)
    #display(ystar)
    x0 = dis*randn(n)
    y0 = dis*randn(m)

    TYPE=""
    A = randn(m,n)
    B = zeros(n,n)
    C = zeros(m,m)
    if prob_type ==1        #quadratic
        println("The Problem is strongly convex-concave quadratic")
        B = random_PD(n)
        C = random_PD(m)
        TYPE="Strongly convex-concave,"

    elseif prob_type ==0   #bilinear
        println("The Problem is bilinear")
        TYPE = "Bilinear,"
    else   # ill-condined quadratic
        B = random_PSD_cond(n, rec_cond)
        C = random_PSD_cond(m, rec_cond)
        A = zeros(m,n)#random_rectangle_cond(n,m, rec_cond)
        TYPE="Ill-cond. quadratic convex-concave,"
    end

    sp  = Saddle_Point(B,A,C,xstar,ystar)
    obj = saddle_point_objective(sp)

    if prob_type ==1
        nabla_F  = obj.∇F
        #nabF_eig = eigen(nabla_F)
        #@show norm(nabF_eig.values[1])
        #@show norm(nabF_eig.values[end])
        #@assert norm(L_eig.values[end]) < 1E-10
        valsB,vecsB = eigen(B)
        @show B_smallest_eig = valsB[1]
        @show B_largest_eig = valsB[end]
        valsC,vecsC = eigen(C)
        @show C_smallest_eig = valsC[1]
        @show C_largest_eig = valsC[end]
    end
