# Some of the solvers in this file are assuming that the Hessian (Jacobian) is constant. See sp_function.jl for the function definiton
"Stepsize smaller than 1e-8 is too small"
stepsize_tol = 1e-8

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
    #ngz = (LinearAlgebra.norm(obj.∇xL(x,y))^2 + LinearAlgebra.norm(obj.∇yL(x,y)^2)^0.5
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
function Armijo_ls(x,y,obj, F, normF, dir, c1)
    beta = 0.5
    inner_max_it = 30
    t = 1
    dir_x = dir[1:n]
    dir_y = dir[n+1:n+m]
    x_p = x +t*dir_x
    y_p = y +t*dir_y
    F_p = [obj.∇xL(x_p,y_p); -obj.∇yL(x_p,y_p)]
    normF_p = LinearAlgebra.norm(F_p)
    "NOTE: I picked the RHS in the 1st condition aribtrarily!"
    RHS = c1*normF
    #println("Right hand side is =$RHS")
    in_it = 0

    while  normF -  normF_p < RHS && in_it < inner_max_it
        LHS = normF -  normF_p
        #println("The left hand side  =$LHS")
        t = beta*t
        x_p = x +t*dir_x
        y_p = y +t*dir_y
        F_p = [obj.∇xL(x_p,y_p); -obj.∇yL(x_p,y_p)]
        normF_p = LinearAlgebra.norm(F_p)
        in_it = in_it + 1
    end
    return t
end
#######################################
function EGM(x,y,obj,dummy0,sp,max_it,prt, inpt, use_ls)
    nablaF = [sp.B  sp.A'
            -sp.A    sp.C]
    eta  = 1/LinearAlgebra.norm(nablaF)
    println("############ Extera Gradient Method ###############")
    println("EGM is using inverse of the Lipschitz constant for stepsize.")
    if inpt==1
        println("Renewing the initial point")
        x = randn(n)
        y = randn(m)
    end
    val, ngx, ngy =0,0,0
    iter=0
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    while ( normF> tol ) && iter < max_it
        xk = x -eta*obj.∇xL(x,y)
        yk = y +eta*obj.∇yL(x,y)

        x= x -eta*obj.∇xL(xk,yk)
        y= y +eta*obj.∇yL(xk,yk)

        F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
        append!(normFAll, normF)
        normF = LinearAlgebra.norm(F)

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
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

"Looks like the bad updating is not working for us. So I omitied from the below function"
function Broyden(x,y,obj,gam,sp,itNum, prt,inpt, use_ls)# gamma is the fixed stepsize
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
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    while ( normF> tol ) && iter < itNum
        z = [x;y]
        p = -H*F
        if(use_ls ==1)
            gam = Armijo_ls(x,y,obj,F,normF,p,c1)
            if gam >= stepsize_tol
                if prt==1
                    println("t = $gam")
                end
                s = gam*p
                z = z+s
                x = z[1:n]
                y = z[n+1:n+m]
                F_old = F
                F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
                Y = F - F_old # Y is y in my notes
                append!(normFAll, normF)
                normF = LinearAlgebra.norm(F)
             else # null step
                if prt==1
                    println("null step")
                end
                s = p
                z_temp = z+s
                x_temp = z_temp[1:n]
                y_temp = z_temp[n+1:n+m]
                F_temp = [obj.∇xL(x_temp,y_temp); -obj.∇yL(x_temp,y_temp)]
                Y = F_temp - F_old
                append!(normFAll, normF)
            end
        end

        stH = s'*H
        H = H + ((s-H*Y)*stH)/(stH*Y)
        #H = H + ((s-H*Y)*Y')/(Y'*Y) %bad update
        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        nx = LinearAlgebra.norm(x)
        ny = LinearAlgebra.norm(y)
        normH = LinearAlgebra.norm(H)

        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("|x| = $nx")
            println("|y| = $ny")
            println("|H^{-1}| = $normH")
            println("iter $iter")
        end
    end
    return x,y,iter,normFAll, val, ngx, ngy

end
#######################################
function secantUpdate_alg(x,y,obj,dummy,sp, itNum,prt,inpt, use_ls)# gamma is the fixed stepsize
    val, ngx, ngy =0,0,0
    println("################ Using Secant Updating with M = I")
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
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    while ( normF> tol ) && iter < itNum
        #println("@@@@@@@@@@@@@@ iter begins")
        z = [x;y]
        H = pinv(D)
        p = -H*F

        #display(D)
        #display(H)
        if(use_ls ==1)
            gama = Armijo_ls(x,y,obj,F,normF,p,c1)
            if gama >= stepsize_tol
                if prt==1
                    println("t = $gama")
                end
                s = gama*p
                z = z+s
                x = z[1:n]
                y = z[n+1:n+m]
                F_old = F
                F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
                Y = F - F_old # Y is y in my notes
                append!(normFAll, normF)
                normF = LinearAlgebra.norm(F)
                yts = Y'*s
                println("y's = $yts")
             else # null step
                 if prt==1
                     println("null step")
                 end
                s = p
                z_temp = z+s
                x_temp = z_temp[1:n]
                y_temp = z_temp[n+1:n+m]
                F_temp = [obj.∇xL(x_temp,y_temp); -obj.∇yL(x_temp,y_temp)]
                Y = F_temp - F_old
                append!(normFAll, normF)
                yts = Y'*s
                println("y's = $yts")
            end
        else #lineseach is not being used
            s = p
            z = z+s
            x = z[1:n]
            y = z[n+1:n+m]
            F_old = F
            F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
            Y = F - F_old # Y is y in my notes
            append!(normFAll, normF)
            normF = LinearAlgebra.norm(F)
            yts = Y'*s
            println("y's = $yts")
        end

        r = Y - D*s # here Ds is just gamma*F_old

        norms2 = LinearAlgebra.norm(s)^2
        "The following two expresion are identical, one is using r the other is using  Y-Ds"
        #D = D+(r*s'+(J*s)*(r'*J)-((s'*J*r)*(J*s*s'))/norms2)/norms2
        E = (Y*s' + J*s*Y'*J-D*s*s' -J*s*s'*D'*J -(s'*J*Y-s'*J*D*s)*J*s*s'/norms2)/norms2
        #rankE = rank(E, rtol=1e-10)
        #println("Secant update rank = $rankE")
        D = D+ E

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
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
function secantShermanWoodbury_alg2(x,y,obj,dummy,sp, itNum,prt,inpt, use_ls)# gamma is the fixed stepsize
    val, ngx, ngy =0,0,0
    println("################ Using Secant Updating with Sherman Woodbury Identity with M = I")
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
    #D = eye
    H = eye
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    while ( normF> tol ) && iter < itNum
        z = [x;y]
        p = -H*F

        if(use_ls ==1)
            gama = Armijo_ls(x,y,obj,F,normF,p,c1)
            if gama >= stepsize_tol
                if prt==1
                    println("t = $gama")
                end
                s = gama*p
                z = z+s
                x = z[1:n]
                y = z[n+1:n+m]
                F_old = F
                F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
                Y = F - F_old # Y is y in my notes
                append!(normFAll, normF)
                normF = LinearAlgebra.norm(F)
                yts = Y'*s
                println("y's = $yts")
             else # null step
                if prt==1
                    println("null step")
                end
                gama=1
                F_old = F
                s = p
                z_temp = z+s
                x_temp = z_temp[1:n]
                y_temp = z_temp[n+1:n+m]
                F_temp = [obj.∇xL(x_temp,y_temp); -obj.∇yL(x_temp,y_temp)]
                Y = F_temp - F_old
                append!(normFAll, normF)
                yts = Y'*s
                println("y's = $yts")
            end
        else #lineseach is not being used
            gama=1
            s = p
            z = z+s
            x = z[1:n]
            y = z[n+1:n+m]
            F_old = F
            F = [obj.∇xL(x,y); -obj.∇yL(x,y)]

            Y = F - F_old # Y is y in my notes
            append!(normFAll, normF)
            normF = LinearAlgebra.norm(F)
            yts = Y'*s
            println("y's = $yts")
        end

        #r = Y - D*s # here Ds is just gamma*F_old
        "since Ds=-stepsize*F_old we get"

        r = Y +gama*F_old

        ns2 = LinearAlgebra.norm(s)^2
        Js = J*s
        Jr = J*r
        α = s'*Jr
        a = r - α*Js/ns2
        Ha = H*a
        Ainv = H - (Ha*s'*H)/(ns2 +s'*Ha )
        AinvJs = Ainv*Js
        H = Ainv - (AinvJs*Jr'*Ainv)/(ns2 + Jr'*AinvJs)

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
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
# NOT calling these functions.
#######################################
function H_secant_alg(x,y,obj,gama,sp, itNum,prt,inpt, use_ls)# gamma is the fixed stepsize
    val, ngx, ngy =0,0,0
    println("################ Using H Secant Updating with M = I")
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
    H = eye

    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    while ( normF> tol ) && iter < itNum
        #println("@@@@@@@@@@@@@@ iter begins")
        z = [x;y]
        #H = pinv(D)
        p = -H*F
        #display(D)
        #display(H)
        if(use_ls ==1)
            gama = Armijo_ls(x,y,obj,F,normF,p,c1)
            if gama >= stepsize_tol
                if prt==1
                    println("t = $gama")
                end
                s = gama*p
                z = z+s
                x = z[1:n]
                y = z[n+1:n+m]
                F_old = F
                F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
                Y = F - F_old # Y is y in my notes
                append!(normFAll, normF)
                normF = LinearAlgebra.norm(F)
             else # null step
                 if prt==1
                     println("null step")
                 end
                s = p
                z_temp = z+s
                x_temp = z_temp[1:n]
                y_temp = z_temp[n+1:n+m]
                F_temp = [obj.∇xL(x_temp,y_temp); -obj.∇yL(x_temp,y_temp)]
                Y = F_temp - F_old
                append!(normFAll, normF)
            end
        end

        #r = s - H*y # here Ds is just gamma*F_old
        normy2 = LinearAlgebra.norm(Y)^2
        E = (s*Y' + J*Y*s'*J-H*Y*Y' -J*Y*Y'*H'*J -(Y'*J*s-Y'*J*H*Y)*J*Y*Y'/normy2)/normy2
        #rankE = rank(E, rtol=1e-10)
        #println("Secant update rank = $rankE")
        H = H+ E

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
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
    normF = LinearAlgebra.norm(F)
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


###############################Line Search methods#############
###############################################################
