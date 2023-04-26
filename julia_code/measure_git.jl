# Measure diversity, difference, and width at sample regions
using DynamicGrids, Crayons, StatsBase, Plots, DataFrames, GLM, Random, JLD
Random.seed!(1234)
nrep=10
size = 129
# Load simulated result of 10 copies 129X129 system
landrept = load("LH2D129-4.jld")["data"]

# Mean diversity/width over all i*i samples in certain copy
smeandiv= Array{Float64,2}(undef,size-1,nrep)
smeanwidth = Array{Float64,2}(undef,size-1,nrep)
# Mean diversity/width over all copies
repmeandiv= Array{Float64,1}(undef,size-1)
repmeanwidth = Array{Float64,1}(undef,size-1)
# Diversity/width at sample regions
diversity = Array{Any,3}(missing,size-1,size-1,size-1)
width = Array{Any,3}(missing,size-1,size-1,size-1)

# Mean min/max over all i*i sample regions in certain copy
smeanlow= Array{Float64,2}(undef,size-1,nrep)
smeanhigh = Array{Float64,2}(undef,size-1,nrep)
# Mean min/max over all copies
repmeanlow= Array{Float64,1}(undef,size-1)
repmeanhigh = Array{Float64,1}(undef,size-1)
# Min/max at sample regions
low = Array{Any,3}(missing,size-1,size-1,size-1)
high = Array{Any,3}(missing,size-1,size-1,size-1)



@time begin
# Iterate over all possible sample size
Threads.@threads for i in (1:size-1)
    println(i)
    # Iterate over 10 copies
    for n in 1:nrep
        # Effective population landscape of certain copy
        landscape=log.(landrept[:,:,n])
    # Take samples at every i X i regions.
    for j1 in 1:size-i, j2 in 1:size-i
        srange1 = (j1:j1+i)
        srange2 = (j2:j2+i)
        samp = @view(landscape[srange1,srange2])
        # Measure diversity/width/min/max of the sample region, store as the (j1,j2,i) element of diversity/width/low/high matrix
        diversity[j1,j2,i]=length(unique(samp))-1
        width[j1,j2,i] = var(samp)*(i^2+2*i)/(i+1)^2
        low[j1,j2,i]=findmin(samp)[1]
        high[j1,j2,i] = findmax(samp)[1]
end
# Average over all size i samples
smeandiv[i,n]=mean(skipmissing(@view(diversity[:,:,i])))
smeanwidth[i,n]=mean(skipmissing(@view(width[:,:,i])))
smeanlow[i,n]=mean(skipmissing(@view(low[:,:,i])))
smeanhigh[i,n]=mean(skipmissing(@view(high[:,:,i])))
end
# Average over all copies of simulations
repmeandiv[i]= mean(@view(smeandiv[i,:]))
repmeanwidth[i]= mean(@view(smeanwidth[i,:]))
repmeanlow[i]= mean(@view(smeanlow[i,:]))
repmeanhigh[i]= mean(@view(smeanhigh[i,:]))
end
end
repmeandiff = repmeanhigh.-repmeanlow;

save("LH2Drepmeandiv129-4.jld","data",repmeandiv)
save("LH2Drepmeanwidth129-4.jld","data",repmeanwidth)
save("LH2Drepmeandiff129-4.jld","data",repmeandiff)
