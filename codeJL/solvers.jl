# The solvers in this file are assuming that the Hessian (Jacobian) is constant. See sp_function.jl for the function definiton
function EGM(x,y,obj,eta,dummy,dummy2,max_it,prt, inpt)
    println("\n############ Extera Gradient Method ###############")
    if inpt==1
        println("Renewing the initial point")
        x = randn(n)
        y = randn(m)
    end
    val, ngx, ngy =0,0,0
    iter=0
    normF=LinearAlgebra.norm(obj.∇xL(x,y)) + LinearAlgebra.norm(obj.∇yL(x,y))
    normFAll=[]
    while ( normF> tol ) && iter < max_it
        xk = x -eta*obj.∇xL(x,y)
        yk = y +eta*obj.∇yL(x,y)

        x= x -eta*obj.∇xL(xk,yk)
        y= y +eta*obj.∇yL(xk,yk)

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        append!(normFAll, normF)
        normF = ngx + ngy
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("iter $iter")
        end
    end
    return x,y,iter,normFAll, val, ngx, ngy
end

#######################################

function Newton_method(x,y,obj,sp)
    println("\n############ Newton Method ###############")
    Hess = [sp.B  sp.A'
            sp.A    -sp.C]
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
    while LinearAlgebra.norm(gz) > tol
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

#######################################
function Resolvent_alg(x,y,obj,gamma,sp)# gamma is the fixed stepsize
    println("\n ############ Using Proximal Point ############")
    nablaF = [sp.B  sp.A'
            -sp.A    sp.C]

    eye = Matrix{Float64}(I, size(nablaF))
    iter=0
    ngz = LinearAlgebra.norm(obj.∇xL(x,y)) + LinearAlgebra.norm(obj.∇yL(x,y))
    while (ngz > tol ) && iter < 5*1e2
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
"Looks like the bad updating is not working for us. So I omitied from the below function"
function Broyden(x,y,obj,gam,sp, M_opt,itNum, prt,inpt)# gamma is the fixed stepsize
    gam = 1
    val, ngx, ngy =0,0,0
    println("############### Broyden's method (good)#####################")
    if inpt==1
        println("Renewing the initial point")
        x = randn(n)
        y = randn(m)
    end
    eye = Matrix{Float64}(I, (n+m,n+m))
    iter=0
    H = eye

    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    normF=LinearAlgebra.norm(obj.∇xL(x,y)) + LinearAlgebra.norm(obj.∇yL(x,y))
    normFAll=[]
    while ( normF> tol ) && iter < itNum
        z = [x;y]
        p = -H*F
        s = gam*p
        z = z+s
        x = z[1:n]
        y = z[n+1:n+m]
        F_old = F
        F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
        Y = F - F_old # Y is y in my notes
        # if M_opt==1
        stH = s'*H
        H = H + ((s-H*Y)*stH)/(stH*Y)
        # else
        #     H = H + ((s-H*Y)*Y')/(Y'*Y)
        # end
        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        append!(normFAll, normF)
        normF = ngx + ngy
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("iter $iter")
        end
    end
    return x,y,iter,normFAll, val, ngx, ngy

end


#######################################
function secantUpdate_alg(x,y,obj,gama,sp, Mopt, itNum,prt,inpt)# gamma is the fixed stepsize
    gama = 1
    val, ngx, ngy =0,0,0
    if Mopt==1
        println("################ Using Secant Updating with M = I")
    elseif Mopt==2
        println("################ Using Secant Updating with M = 0.5(D+D')")
    else
        println("################ Using Secant Updating with M = 0.5(D_+ + D_+')")
    end
    if inpt==1
        println("Renewing the initial point")
        x = randn(n)
        y = randn(m)
    end

    eye1 = Matrix{Float64}(I, size(sp.B))
    eye2 = Matrix{Float64}(I, size(sp.C))
    zeroBlock = zeros(size(sp.A))
    sizezeroBlock = size(sp.A)
    J = [eye1       zeroBlock'
        zeroBlock      -eye2]
    eye = Matrix{Float64}(I, size(J))
    iter=0
    D = eye

    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    normF=LinearAlgebra.norm(obj.∇xL(x,y)) + LinearAlgebra.norm(obj.∇yL(x,y))
    normFAll=[]
    while ( normF> tol ) && iter < itNum
        z = [x;y]
        #normD = LinearAlgebra.norm(D)
        H = pinv(D)
        #normH = LinearAlgebra.norm(H)
        #println("normD = $normD")
        #println("normH = $normH")
        p = -H*F
        s = gama*p
        z = z+s
        x = z[1:n]
        y = z[n+1:n+m]
        F_old = F
        F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
        Y = F - F_old # Y is y in my notes
        r = Y - D*s # here Ds is just gamma*F_old
        if Mopt==1
            " Assuming M = I:"
            norms2 = LinearAlgebra.norm(s)^2
            "The following two expresion are identical, one is using r the other is using  Y-Ds"
            #D = D+(r*s'+(J*s)*(r'*J)-((s'*J*r)*(J*s*s'))/norms2)/norms2
            D = D+ (Y*s' + J*s*Y'*J-D*s*s' -J*s*s'*D'*J -(s'*J*Y-s'*J*D*s)*J*s*s'/norms2)/norms2
        elseif Mopt==2
            "Assuming M = 0.5(D+D'):"
            M = 0.5*(D+D')
            norms2 = (s'*M*s)
            D = D+(r*s'*M +M*(J*s)*(r'*J)-((s'*J*r)*(M*(J*s)*(s'*M)))/norms2)/norms2
            #D = D + (Y*(s'*M)+M*(J*s)*(Y'*J) -D*s*s'*M -M*J*s*s'*D'*J -(s'*J*Y -s'*J*D*s)*(M*J*s*s'*M)/norms2)/norms2

        else
            " Assuming M = 0.5(D_+ + D_+'):"
        end

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        append!(normFAll, normF)
        normF = ngx + ngy
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("iter $iter")
        end
    end
    return x,y,iter,normFAll, val, ngx, ngy
end

#######################################
function secantUpdate_extra_grad_alg(x,y,obj,gamma,sp, itNum)# gamma is the fixed stepsize
    println("\n ############ Using Secant Update With Extra Gradient ############")
    val, ngx, ngy =0,0,0
    eye1 = Matrix{Float64}(I, size(sp.B))
    eye2 = Matrix{Float64}(I, size(sp.C))
    zeroBlock = zeros(size(sp.A))
    sizezeroBlock = size(sp.A)
    J = [eye1       zeroBlock'
        zeroBlock      -eye2]
    eye = Matrix{Float64}(I, size(J))
    iter=0
    D = eye

    xk = x -gamma*obj.∇xL(x,y)
    yk = y +gamma*obj.∇yL(x,y)
    Fk = [obj.∇xL(xk,yk) ; -obj.∇yL(xk,yk)]
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    normF = LinearAlgebra.norm(obj.∇xL(x,y)) + LinearAlgebra.norm(obj.∇yL(x,y))
    while (normF > tol ) && iter < itNum
        z = [x;y]
        #Instead using gradeint of (x,y) use (xk,yk)
        p = -pinv(D)*Fk
        s = gamma*p
        z = z+s
        x = z[1:n]
        y = z[n+1:n+m]
        F_old = F
        xk = x -gamma*obj.∇xL(x,y)
        yk = y +gamma*obj.∇yL(x,y)
        Fk = [obj.∇xL(xk,yk) ; -obj.∇yL(xk,yk)]
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
    return x,y,iter, normFAll, val, ngx, ngy
end
