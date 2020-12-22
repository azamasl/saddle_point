struct Saddle_Point
    # L = 0.5(x-xstar)'*B*(x-xstar) +(y-ystar)'*A*(x-xstar) -0.5(y-ystar)'*C*(y-ystar)
    B::Array{Float64,2} # PSD matrix
    A::Array{Float64,2}
    C::Array{Float64,2} # PSD matrix
    xstar::Array{Float64,1}
    ystar::Array{Float64,1}
end

struct ObjectiveFunction
    L::Function # objective function
    ∇xL::Function # (sub)gradient of objective
    ∇yL::Function # (sub)gradient of objective
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

    return ObjectiveFunction(L,∇xL,∇yL)
end

function sad_point2D()
    function L(x,y)
        return  (x^2-1)*(x^2-9)+x*y -(y^2-1)*(y^2-9)
    end
    function ∇xL(x,y)
        return 4*x^3-20*x + y
    end
    function ∇yL(x,y)
        return 4*y^3-20*y + x
    end
    return ObjectiveFunction(L,∇xL,∇yL)
end



####################### Generating random matrices
"Generate random matrix ∈ S^n with entries ~ N(0,1)"
function random_sym(N)
    S = randn(N, N)
    return (S + S')/sqrt(2)
end
"Generate random PSD matrix with entries: off-diagonal entries ∼ N(0,1), diagonal entries ~ N(√n,2)"
function random_PSD(N)
    S = randn(N, N) / (N ^ 0.25)
    S = S * S'
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
