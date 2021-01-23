using  Plots, JLD2
include("problem_settings.jl")
ymin, ymax = 1E-10,1E20

@load("output/ws_tr_0 _ 220.jld2")

nfs = ret.nfs
plot(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[1])

@load("output/ws_EGM_0 _ 220.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[2])

@load("output/ws_lsSec_0 _ 220.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[4])

plot!(xlabel = "Iteration", ylabel = "||F||", title = string(TYPE," m = $m, n=$n")
,legend = :topleft,legendfont = font(7), ylims=(ymin,ymax))

savefig("plots/all_$prob_type.png")
