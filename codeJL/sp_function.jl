struct Saddle_Point
    # L = (x-xstar)'*B*(x-xstar) +(y-ystar)'*A*(x-xstar) -(y-ystar)'*C*(y-ystar)
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
        return  (x-xstar)'*B*(x-xstar) +(y-ystar)'*A*(x-xstar) -(y-ystar)'*C*(y-ystar)
    end
    function ∇xL(x,y)
        return 2*B*(x-xstar) + A'*(y-ystar)
    end
    function ∇yL(x,y)
        return A*(x-xstar) - 2*C'*(y-ystar)
    end
    return ObjectiveFunction(L,∇xL,∇yL)
end

####################### Generating random matrices
"Generate random matrix ∈ S^n with entries ~ N(0,1)"
function random_sym(n)
    S = randn(n, n)
    return (S + S')/sqrt(2)
end
"Generate random PSD matrix with entries: off-diagonal entries ∼ N(0,1), diagonal entries ~ N(√n,2)"
function random_PSD(n)
    S = randn(n, n) / (n ^ 0.25)
    S = S * S'
end
