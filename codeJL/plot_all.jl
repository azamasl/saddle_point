using  Plots, JLD2
ymin, ymax = 1E-10,1E10


lb = ["Tr-Secant" "EGM" "Broyden" "Ls-Secant" "NoLs-Secant"]
#lb = ["Tr-Secant"  "Ls-Secant" "NW41-Secant"]

@load("output/ws_tr_0-900.jld2")
nfs = ret.nfs
plot(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[1])


@load("output/ws_EGM_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[2])


@load("output/ws_BROYDEN_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[3])


@load("output/ws_lsSec_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[4])


@load("output/ws_noSec_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[5])




# @load("output/ws_tr_NW41_0-900.jld2")
# nfs = ret.nfs
# plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[3])







plot!(xlabel = "Iteration", ylabel = "||F||", title = string(TYPE," m = $m, n=$n")
,legend = :topright,legendfont = font(7), ylims=(ymin,ymax))

savefig("plots/all_$prob_type.png")
