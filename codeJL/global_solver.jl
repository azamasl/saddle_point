"Computes the second direction length in the dogleg dirc."
function getAlpha(p,q,Del)
    "p'p + alpha^2q'q + 2alpha p'q = Del^2  "
    a = q'*q
    b = 2*p'*q
    c = p'*p-Del^2
    #TODO: make sure the alternative alpha is checked too
    alpha = (-b + sqrt(b^2 -4*a*c))/(2*a)
    alpha2 = (-b - sqrt(b^2 -4*a*c))/(2*a)
    print("b = $b ")
    println("roots (second half of the dogled direction) are : $alpha and $alpha2 ")
    return alpha
end


"Uses dogleg dir in a trust-region framework"
function tr_dogleg(x,y, obj,sp, itNum,prt, F_tol, Del, max_Del, eta)
    ################First let ∇F be constant:
    n = length(x)
    m = length(y)

    nabla_F = I(m+n)
    if (m==1 && n==1)
        nabla_F  = obj.∇F(x,y)
     else
    #     nabla_F = [sp.B    sp.A'
    #               -sp.A    sp.C]
        nabla_F  = obj.∇F
    end

    #eta  = 1/LinearAlgebra.norm(nablaF)
    #@show nabla_F
    val, ngx, ngy =0,0,0
    println("########### Inside the tr_dogleg method :")
    #@show cond(sp.B)
    #@show cond(sp.A)
    #@show cond(sp.C)
    @show cond(nabla_F)
    J = [I(n)       zeros(n,m)
        zeros(m,n)     -I(m)]
    B_k = I(n+m)
    H_k = I(n+m)

    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    normF=norm(F)
    g = nabla_F'*F
    k=0    # count all iterations
    it=0   # count only actual steps

    normFAll=[]
    while ( normF> F_tol ) && k < itNum
        z = [x;y]
        "compute the direction"
        pB = -H_k*(H_k'*g)                 #"compute quasi-Newton step"
        if( norm(pB) <=Del)
            s_k= pB
        else                               #compute the Cauchy point
            Bg = B_k*g
            norm_g = norm(g)
            Del2n_g=Del/norm_g
            ps = -Del2n_g*g
            norm_Bg_sqr = Bg'*Bg
            n_g2n_Bg=(norm_g)^2/norm_Bg_sqr
            #tau = n_g2n_Bg/Del2n_g
            #if tau >= 1                  # Cauchy point outside
            if n_g2n_Bg >= Del2n_g        # Cauchy point outside
                s_k = ps
            else                          # Cauchy point inside, move along dogleg
                pU = -n_g2n_Bg*g
                q = pB-pU
                alpha = getAlpha(pU,q,Del)
                s_k = pU + alpha*q
            end
        end
        "compue the prediction"
        z_new = z+s_k
        F_new = [obj.∇xL(z_new[1:n],z_new[n+1:n+m]); -obj.∇yL(z_new[1:n],z_new[n+1:n+m])]
        numer = 0.5*(norm(F)^2 - norm(F_new)^2)
        Bs = B_k*s_k
        denom = -g'*s_k - 0.5*(Bs)'*(Bs)
        if denom ==0
            rho = 1e99                 #division by 0
        else
            rho = numer/denom
        end
        "expand or contract"
        if rho < 0.25                 #contract
            Del = 0.5*Del
        else                          #expand
            Del = min(2*Del, max_Del)
        end
        "update B and H"
        Y = F_new - F
        r = Y - Bs
        norms2 = norm(s_k)^2
        #β = 2*rand()                  # we prefer β=1 whenever B_{k+1} is nonsingular, which is almost surely the case
        β = 1
        if prt==1
            println("β = $β")
        end
        B_k = B_k + β*(J*s_k*r'*J+ r*s_k')/norms2  - β^2*(s_k'*J*r*J*s_k*s_k')/(norms2)^2
        H_k = inv(B_k)

        if rho >=eta # sufficient decrease
            z = z_new
            x = z_new[1:n]
            y = z_new[n+1:n+m]
            F = F_new
            append!(normFAll, normF)
            normF = LinearAlgebra.norm(F)
            if (m==1 && n==1)
                nabla_F  = obj.∇F(x,y)
            end
            g = nabla_F'*F

            val = norm(obj.L(x,y))
            ngx = norm(obj.∇xL(x,y))
            ngy = norm(obj.∇yL(x,y))
            it = it+1
        else
            if prt==1
                println("the step was null")
            end
        end
        k=k+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            ng = sqrt(ngx^2+ngy^2)
            println("||g|| = $ng")
            println("#################################End of iter $k")
        end
    end#while loop
    ng = sqrt(ngx^2+ngy^2)
    return x,y,it,normFAll, val, ng,k
end
