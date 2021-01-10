using NLsolve
F(x) = [4*x[1]^3-20*x[1]+x[2]
          4*x[2]^3-20*x[2]-x[1]]

#results = nlsolve(F, [ rand()-rand(); rand()-rand()])
#for i=1:9
    @show print(nlsolve(F, [2*(rand()-rand()); 2*(rand()-rand()) ]))
    #contour(F)
#end

"computing the zeros analytically:"
#x = 4*y*(y^2-5)
#y = 4*(4*y*(y^2-5))^3 -80*y*(y^2-5)+y = 0

using Roots
f(y) = 4*(4*y*(y^2-5))^3 -80*y*(y^2-5)+y
#@show print(fzero(f, 0,1000))
@show f(-0.112226109193021)
