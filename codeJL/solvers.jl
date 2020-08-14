# The solvers in this file are assuming that the Hessian (Jacobian) is constant. See sp_function.jl for the function definiton
function EGM_solver(x,y,obj,eta)
    iter=0
    while LinearAlgebra.norm(obj.∇xL(x,y)) > 1e-6 || LinearAlgebra.norm(obj.∇yL(x,y)) > 1e-6
        xk = x -eta*obj.∇xL(x,y)
        yk = y +eta*obj.∇yL(x,y)

        x= x -eta*obj.∇xL(xk,yk)
        y= y +eta*obj.∇yL(xk,yk)
        #val =obj.L(x,y)
        #ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        #ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        #println("L = $val")
        #println("|∇xL| = $ngx")
        #println("|∇yL| = $ngy")
        iter=iter+1
    end
    return x,y,iter
end
#######################################
function Resolvent_alg(x,y,obj,gamma,sp)# gamma is the fixed stepsize
    nablaF = [2*sp.B  sp.A'
            -sp.A    2*sp.C]

    eye = Matrix{Float64}(I, size(nablaF))
    iter=0
    while (LinearAlgebra.norm(obj.∇xL(x,y)) > 1e-8 || LinearAlgebra.norm(obj.∇yL(x,y)) > 1e-8 ) && iter < 5*1e2
        z = [x;y]
        F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
        z = z-gamma*pinv(eye+gamma*nablaF)*F
        x = z[1:n]
        y = z[n+1:n+m]
        # val =obj.L(x,y)
        # ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        # ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        # println("L = $val")
        # println("|∇xL| = $ngx")
        # println("|∇yL| = $ngy")
        iter=iter+1
    end
    return x,y,iter
end

#######################################
function secantUpdate_alg(x,y,obj,gamma,sp)# gamma is the fixed stepsize
    eye1 = Matrix{Float64}(I, size(sp.B))
    eye2 = Matrix{Float64}(I, size(sp.C))
    zeroBlock = zeros(size(sp.A))
    sizezeroBlock = size(sp.A)
    J = [eye1       zeroBlock'
        zeroBlock      -eye2]
    eye = Matrix{Float64}(I, size(J))
    iter=0
    D = eye
    #F_old = [ zeros(size(obj.∇xL(x,y))); zeros(size(obj.∇yL(x,y)))]
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    while (LinearAlgebra.norm(obj.∇xL(x,y)) > 1e-8 || LinearAlgebra.norm(obj.∇yL(x,y)) > 1e-8 ) && iter < 2*1e3
        z = [x;y]
        p = -pinv(D)*F
        s = gamma*p
        z = z+s
        x = z[1:n]
        y = z[n+1:n+m]
        F_old = F
        F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
        Fdiff = F - F_old # Fdiff is y in my notes
        r = Fdiff - D*s # here Ds is just gamma*F_old
        norms2 = LinearAlgebra.norm(s)^2
        D = D+(r*s'+(J*s)*(r'*J)-((s'*J*r)*(J*s*s'))/norms2)/norms2
        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        println("L = $val")
        println("|∇xL| = $ngx")
        println("|∇yL| = $ngy")
        iter=iter+1
        println("end of iter $iter")
    end
    return x,y,iter
end


#######################################

function Newton_method(x,y,obj,sp)
    Hess = [2*sp.B  sp.A'
            sp.A    -2*sp.C]
    #Q,R = qr(Hess)
    #if min(abs(diag(R))) == 0 #cols of Hess are independent iff diag(R) has no 0
    invHess = pinv(Hess) #when Hess is invertible pinv is just inv
    #invHess = inv(Hess)
    #eye = Matrix{Float64}(I, size(Hess))
    #invHess = Hess\eye

    iter=0
    z = [x;y]
    gz = [obj.∇xL(x,y);obj.∇yL(x,y)]
    n = size(x,1)
    m = size(y,1)
    while LinearAlgebra.norm(gz) > 1e-8
    #for k=1:10
        z = z -invHess*gz
        x = z[1:n]
        y = z[n+1:n+m]
        gz = [obj.∇xL(x,y);obj.∇yL(x,y)]
        #val =obj.L(x,y)
        #ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        #ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        # println("L = $val")
        # println("|∇xL| = $ngx")
        # println("|∇yL| = $ngy")
        iter=iter+1
    end
    return x,y,iter
end
