struct Saddle_Point
    # L = 0.5(x-xstar)'*B*(x-xstar) +(y-ystar)'*A*(x-xstar) -0.5(y-ystar)'*C*(y-ystar)
    B::Array{Float64,2} # PSD matrix
    A::Array{Float64,2}
    C::Array{Float64,2} # PSD matrix
    xstar::Array{Float64,1}
    ystar::Array{Float64,1}
    #∇F::Array{Float64,2} # PSD matrix
end

struct fun_return
    x_sol::Array{Float64}
    y_sol::Array{Float64}
    iter::Int
    nfs::Array{Float64}
    val::Float64
    ng::Float64
end


struct ObjectiveFunction
    L::Function # objective function
    ∇xL::Function # (sub)gradient of objective
    ∇yL::Function # (sub)gradient of objective
    ∇F::Array{Float64,2}
end

struct Obje_Fun2D
    L::Function # objective function
    ∇xL::Function # (sub)gradient of objective
    ∇yL::Function # (sub)gradient of objective
    ∇F::Function
end


"saddle-point function L(x,y)= x'Bx+y'Ax-y'Cy with (sub)gradient ∇L(x,y)."
function saddle_point_objective(sp::Saddle_Point)
    B, A,C,xstar,ystar = sp.B, sp.A, sp.C,sp.xstar, sp.ystar
    function L(x,y)
        return  0.5*(x-xstar)'*B*(x-xstar) +(y-ystar)'*A*(x-xstar) -0.5*(y-ystar)'*C*(y-ystar)
    end
    function ∇xL(x,y)
        return B*(x-xstar) + A'*(y-ystar)
    end
    function ∇yL(x,y)
        return A*(x-xstar) - C'*(y-ystar)
    end
    ∇F = [sp.B    sp.A'
          -sp.A    sp.C]
    return ObjectiveFunction(L,∇xL,∇yL,∇F)
end

function sad_point2D_objective()
    function L(x,y)
        return  (x.^2 .- 1).*(x.^2 .- 9) .+ x.*y .- (y.^2 .- 1).*(y.^2 .- 9)
    end
    function ∇xL(x,y)
        return 4*x.^3 .- 20*x .+ y
    end
    function ∇yL(x,y)
        return -4*y.^3 .+ 20*y .+ x
    end
    function ∇F(x,y)
        return  [12*x'*x-20    1.
                  -1.    12*y'*y-20]
    end
    return Obje_Fun2D(L,∇xL,∇yL,∇F)
end

################ random numbers form the interval
function rand_interval(l,u)
    return   l + (u-l)*rand()
end
####################### Generating random matrices
"Generate random matrix ∈ S^n with entries ~ N(0,1)"
function random_sym(N)
    S = randn(N, N)
    return (S + S')/sqrt(2)
end

"Generate random PSD matrix with entries: off-diagonal entries ∼ N(0,1), diagonal entries ~ N(√n,2)"
function random_PSD(N, scale)
    S = randn(N, N) / (N ^ 0.25)
    S = scale*(S * S')#scaling S down to avoid overflow for larger matrices
    return S
end


"Generate random PD matrix with entries: off-diagonal entries ∼ N(0,1), diagonal entries ~ N(√n,2)"
function random_PD(N,scale)
    S = randn(N, N) / (N ^ 0.25)
    S = scale*(S * S')#scaling S down to avoid overflow for larger matrices
    eig_vals, eig_vecs = eigen(S)
    #adding 1e-3 to small eigenvalues to make sure S is PD
    eig_vals[findall(eig_vals .< 1E-3)] =  eig_vals[findall(eig_vals .< 1E-3)] .+ 1E-3
    S = eig_vecs*Diagonal(eig_vals)*eig_vecs'
    return S
end

"Generate random matrix with  entries ∼ N(0,1)"
function random_scaled(M,N,scale)
    S = scale*randn(M, N)
    return S
end


"Generates a uniformly random PSD matrix with the given Reciprocal Condition Mumber: rcond"
#$$\sigma_1\Big( 1-       \frac{(1 - \frac{1}{c})}{\sigma_1-\sigma_N}(\sigma_1-s)   \Big)   $$
function random_PSD_cond(N, rcond)
    Mat = 2 .* rand(N,N).-1 #  -1 < entries <1
    #Mat = (Mat+Mat')/2
    U,s,V = svd(Mat) # s is a vector and sorted descendingly
    # shifting&scaling s such that its cond. is c:
    curC = cond(Mat)
    println("current cond number = $curC ")
    #display(s)
    c = 1/rcond
    println("desired cond number =$c ")
    s = s[1] .* (  1 .-   ((1-rcond)/(s[1]-s[end])) .* (s[1] .- s)     )
    Mat = U*Diagonal(s)*U'
    @show cond(Mat)
    ##TODO: remove this, only for sanity check
    #@show isposdef(Mat)
    return Mat
end


#
# "Generates a uniformly random rectangle matrix with the given Reciprocal Condition Mumber: rcond"
# #$$\sigma_1\Big( 1-       \frac{(1 - \frac{1}{c})}{\sigma_1-\sigma_N}(\sigma_1-s)   \Big)   $$
# function random_rectangle_cond(p,q, rcond)
#     Mat = 2 .* rand(p,q).-1 #  -1 < entries <1
#     U,s,Vt = svd(Mat) # s is a vector and sorted descendingly
#     @show size(U)
#     @show size(s)
#     @show size(Vt)
#     # shifting&scaling s such that its cond. is c:
#     curC = cond(Mat)
#     println("current cond number = $curC ")
#     #display(s)
#     c = 1/rcond
#     println("desired cond number =$c ")
#     s = s[1] .* (  1 .-   ((1-rcond)/(s[1]-s[end])) .* (s[1] .- s)     )
#     #@show size(s)
#     Di = zeros(p,q)
#     #ls = length(s)
#     t = min(p,q)
#     for i=1:t
#         Di[i,i] = s[i]
#     end
#     Mat = U*Di*Vt
#     @show cond(Mat)
#     return Mat
# end
