using  Plots, JLD2
ymin, ymax = 1E-10,1E5


lb = ["Tr-Secant" "EGM" "Ls-Broyden" "Ls-Secant" "NoLs-Secant" "NoLs-Broyden"]
#lb = ["Tr-Secant"  "Ls-Secant" "NW41-Secant"]

@load("output/ssc_scale-1/ws_tr_1-900.jld2")
nfs = ret.nfs
plot(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[1])
#

@load("output/ssc_scale-1/ws_EGM_1-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[2])


@load("output/ssc_scale-1/ws_BROYDEN_1-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[3])


@load("output/ssc_scale-1/ws_lsSec_1-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[4])


@load("output/ssc_scale-1/ws_noSec_1-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[5])


@load("output/ssc_scale-1/ws_noBROYDEN_1-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[6])


# @load("output/ssc_scale-1/ws_tr_NW41_1-900.jld2")
# nfs = ret.nfs
# plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[3])







plot!(xlabel = "Iteration", ylabel = "||F||", title = string(TYPE," m = $m, n=$n")
,legend = :bottomright,legendfont = font(7), ylims=(ymin,ymax))

savefig("plots/all_$prob_type.png")
