using LinearAlgebra

function zhat(gamma,nablaF,F,zt)
    eye = Matrix{Float64}(I, size(nablaF))
    return zt-gamma*pinv(eye+gamma*nablaF)*F
end

# Bregman Divergence, assuming it's squared norm 2
function D(z,zt)
    return 0.5*LinearAlgebra.norm(z-zt)^2
end

# Binary Search for stepsize
function bs(z,eps, lo, hi, nablaF,F,T,min_sing,normJ)
    #println("bs is being called")
    mid = (lo+hi)/2
    normF = LinearAlgebra.norm(F)
    delta = min_sing/(12*normF)
    deltainv = 1/delta
    C = deltainv^2*((deltainv+normJ)/(12*(min_sing*normF)))^3
    Cbar = max(C,1)
    N = log(Cbar*T/delta)

    for k=0:N-1
        zh = zhat(mid,nablaF,F,z)
        D = 1/(12*LinearAlgebra.norm(zh-z))
        if mid <= D
            lo = mid
        else
            hi = mid
        end
        mid=(lo+hi)/2
    end
    return hi
end



function alg2(zt,eps,T,fun,sp,First_Order)
    Gamma=0
    sum_zhat = zeros(n+m,1)
    wei_zhat = zeros(n+m,1)
    nablaF =[2*sp.B  sp.A'  #Jacobian of F is constant
            -sp.A    +2*sp.C]
    svals = svdvals(nablaF)#singvales are sorted in descending order
    #l= length(svals)
    min_sing = svals[n+m]
    max_sing = svals[1]
    normJ = LinearAlgebra.norm(nablaF)
    println("min_sing = $min_sing")
    #println("max_sing = $max_sing")
    println("norm = $normJ")
    eye = Matrix{Float64}(I, size(nablaF))
    for t=1:T
        #println("iter $t ")
        x = zt[1:n]
        y = zt[n+1:n+m]
        F = [fun.∇xL(x,y); -fun.∇yL(x,y)]#note the -
        normF = LinearAlgebra.norm(F)
        gamma_lo = min_sing/(12*normF)
        gamma_hi = T^1.5
        zh = zhat(gamma_hi,nablaF,F,zt)
        Dz = 1/(8*LinearAlgebra.norm(zh-zt))
        if gamma_hi < Dz
            gamma = gamma_hi
        elseif gamma_lo >= gamma_hi
            gamma = gamma_lo
        else
            gamma = bs(zt,eps,gamma_lo,gamma_hi, nablaF,F,T,min_sing,normJ)
        end
        zh = zhat(gamma,nablaF,F,zt) # eq. (12)
        #isplay(zh)
        x = zh[1:n]
        y = zh[n+1:n+m]
        F = [fun.∇xL(x,y); -fun.∇yL(x,y)]
        #solving the argmin in (14), assuming D is the squared_norm2

        if First_Order==1 #the reason this update doesn't exactly look like (12) is in (14) they use D(z,zt), instead of D(z,zh)
            zt = zh + pinv(eye+gamma*nablaF)*(zt -zh -gamma*F)
        elseif First_Order==2
            zt = zhat(gamma, nablaF, F, zh) # try with D(z,zh) in first-order (14)
        else # zeros_order, setting the gradient to zero to compute the argmin.
            zt = zt-gamma*F
        end

        Gamma = Gamma+gamma
        sum_zhat = sum_zhat + gamma*zh
        wei_zhat = sum_zhat/Gamma
        #println("Current solution: $wei_zhat")
    end
    return wei_zhat
end
