# Some of the solvers in this file are assuming that the Hessian (Jacobian) is constant. See sp_function.jl for the function definiton
"Stepsize smaller than 1e-8 is too small"
stepsize_tol = 1e-8

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
        #LHS = normF -  normF_p
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
function EGM(x,y,obj,dummy0,sp,max_it,prt, use_ls,Ftol)
    # nablaF = [sp.B  sp.A'
    #         -sp.A    sp.C]
    nablaF  = obj.∇F
    eta  = 1/LinearAlgebra.norm(nablaF)
    println("############ Extra Gradient Method ###############")
    println("EGM is using inverse of the Lipschitz constant for stepsize.")
    val, ngx, ngy =0,0,0
    iter=0
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    while ( normF> Ftol ) && iter < max_it
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
    ng = sqrt(ngx^2+ngy^2)
    return x,y,iter,normFAll, val, ng
end

#######################################
function Broyden(x,y,obj,fixed_stepsz,sp,itNum, prt,use_ls, Ftol)
    println("Broyden's method")
    val, ngx, ngy =0,0,0

    H = I(m+n)#diag_PD#
    if prob_type==0
        println("Borydent method with H0=I is undefined for bilinear problems. In this case we set H0 to a random diagonal PD:")
        vec = randn(n+m)
        vec = vec.^2
        #making sure eich eigenvalue is at least 1E-1.
        vec[findall(vec .< 1E-1)] =  vec[findall(vec .< 1E-1)] .+ 1E-1
        diag_PD = Diagonal(abs.(vec))# random diagonal psd matrix
        #display(diag_PD)
        H = diag_PD
    end
    ########## ONLY for testing purpose
    # D = [sp.B  sp.A'
    #     -sp.A  sp.C]
    # z0= [x;y]
    # #display(Hess)
    # s0ty0 = -z0'*(D*D*D)*z0
    # println("With H_0=I and for bilinear functions, s_0ty_0 must be 0, here it is:  $s0ty0")

    ######################
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    iter=0
    while ( normF> Ftol ) && iter < itNum
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
        else # no line search
            gam=fixed_stepsz
            s = gam*p
            z = z+s
            x = z[1:n]
            y = z[n+1:n+m]
            F_old = F
            F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
            Y = F - F_old # Y is y in my notes
            append!(normFAll, normF)
            normF = LinearAlgebra.norm(F)
        end

        stH = s'*H
        stHy = stH*Y
        sty = s'*Y
        H = H + ((s-H*Y)*stH)/stHy

        val = obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("###End of iter $iter")
        end
    end
    ng = sqrt(ngx^2+ngy^2)
    return x,y,iter,normFAll, val, ng

end

#######################################
#Our secant method.
function secant_inv(x,y,obj,fixed_stepsz,sp, itNum,prt, use_ls,Ftol)# gamma is the fixed stepsize
    val, ngx, ngy =0,0,0    
    println("Line search secant method")
    J = [I(n)       zeros(n,m)
        zeros(m,n)     -I(m)]
    H = I(n+m)

    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    iter=0
    while ( normF> Ftol ) && iter < itNum
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
            end
        else #lineseach is not being used
            gama=fixed_stepsz
            s = gama*p
            z = z+s
            x = z[1:n]
            y = z[n+1:n+m]
            F_old = F
            F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
            Y = F - F_old # Y is y in my notes
            append!(normFAll, normF)
            normF = LinearAlgebra.norm(F)
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
        denom1 =(ns2 +s'*Ha )
        Ainv = H - (Ha*s'*H)/denom1
        AinvJs = Ainv*Js
        denom2 =(ns2 + Jr'*AinvJs)
        H = Ainv - (AinvJs*Jr'*Ainv)/denom2

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("#################################End of iter $iter")
        end
    end
    ng = sqrt(ngx^2+ngy^2)
    return x,y,iter,normFAll, val, ng
end
