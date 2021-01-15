using JLD2

@load("ws0.jld2")

println(" Print result in order for: tr_dog_leg EGM Broyden secant_ls secant_Nols")


for i=1:fun_num:
    println("L = $val[i]")
    println("|∇L| = $ng[i]")
    println("number of iterations  = $iter[i]")
end
