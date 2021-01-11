using  Plots, JLD2

@load("ws0.jld2")
#@load("test.jld2")
l1 = size(nfs[1],1)
l2 = size(nfs[2],1)
l3 = size(nfs[3],1)
l4 = size(nfs[4],1)
l5 = size(nfs[5],1)
ymin, ymax = 1E-10,1E20

plot(range(1,l1,step=1), nfs[1], yscale = :log10, label = lb[1])
plot!(range(1,l2,step=1), nfs[2], yscale = :log10, label = lb[2])
plot!(range(1,l3,step=1), nfs[3], yscale = :log10,label = lb[3])
plot!(range(1,l4,step=1), nfs[4], yscale = :log10,label = lb[4])
plot!(range(1,l5,step=1), nfs[5], yscale = :log10,label = lb[5])
plot!(xlabel = "Iteration", ylabel = "||F||", title = string(TYPE," m = $m, n=$n")
,legend = :bottomleft,legendfont = font(7), ylims=(ymin,ymax))
savefig("sample_plots/globalName.png")
