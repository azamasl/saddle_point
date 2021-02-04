using  Plots, JLD2
ymin, ymax = 1E-10,1E20


lb = ["Tr-Secant" "EGM" "Ls-Broyden" "Ls-Secant" "NoLs-Secant" "NoLs-Broyden"]
#lb = ["Tr-Secant"  "Ls-Secant" "NW41-Secant"]

@load("output/bil_scale-1/ws_tr_0-900.jld2")
nfs = ret.nfs
plot(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[1])
@show tr_iter = ret.iter

@load("output/bil_scale-1/ws_EGM_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[2])
@show EGM_iter = ret.iter

@load("output/bil_scale-1/ws_BROYDEN_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[3])
@show lsBROY_iter = ret.iter

@load("output/bil_scale-1/ws_lsSec_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[4])
@show lsSec_iter = ret.iter

@load("output/bil_scale-1/ws_noSec_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[5])
@show noSec_iter = ret.iter

@load("output/bil_scale-1/ws_noBROYDEN_0-900.jld2")
nfs = ret.nfs
plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[6])
@show noBROY_iter = ret.iter

# @load("output/bil_scale-1/ws_tr_NW41_0-900.jld2")
# nfs = ret.nfs
# plot!(range(1,size(nfs,1),step=1), nfs, yscale = :log10, label = lb[3])







plot!(xlabel = "Iteration", ylabel = "||F||", title = string(TYPE," m = $m, n=$n")
,legend = :bottomright,legendfont = font(7), ylims=(ymin,ymax))

savefig("plots/all_$prob_type.png")
#/Users/azam/google/CurrentWork/minmax/codeJL
