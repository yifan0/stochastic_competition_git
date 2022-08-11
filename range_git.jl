# Measure range at sample regions
using DynamicGrids, Crayons, StatsBase, Plots, DataFrames, GLM, Random, JLD

# 10 independent copies of 129X129 system
size = 129
nrep = 10
# Load simulated result of 10 copies 129X129 system
landrept = load("LH2D129-4.jld")["data"]

# Mean min/max over all i*i sample regions in certain copy
smeanlow= Array{Float64,2}(undef,size-1,nrep)
smeanhigh = Array{Float64,2}(undef,size-1,nrep)
# Mean min/max over all copies
repmeanlow= Array{Float64,1}(undef,size-1)
repmeanhigh = Array{Float64,1}(undef,size-1)

# Min/max at sample regions
low = Array{Any,3}(missing,size-1,size-1,size-1)
high = Array{Any,3}(missing,size-1,size-1,size-1)

# Iterate over all possible sample size
for i in (1:size-1)
    # Iterate over 10 copies
    for n in 1:nrep
        # Effective population landscape of certain copy
        landscape=log.(landrept[:,:,n])
    # Take samples at every s patches. Take (size-i) samples in total.
    s = floor(Int, sqrt(size-i))
    for j1 in 1:s, j2 in 1:s
    # An iXi region (srange1)X(srange2) is taken as a sample
    srange1 = (j1*s+1:j1*s+i)
    srange2 = (j2*s+1:j2*s+i)
    # Measure min/max of the sample region, store as the (j1,j2,i) element of min/max matrix
    low[j1,j2,i]=findmin(landscape[srange1,srange2])[1]
    high[j1,j2,i] = findmax(landscape[srange1,srange2])[1]
end
# Average over all size i samples
smeanlow[i,n]=mean(skipmissing(@view(low[:,:,i])))
smeanhigh[i,n]=mean(skipmissing(@view(high[:,:,i])))
end
# Average over all copies of simulations
repmeanlow[i]= mean(@view(smeanlow[i,:]))
repmeanhigh[i]= mean(@view(smeanhigh[i,:]))
end

# Calculate the range
repmeandiff = repmeanhigh.-repmeanlow;

# Save data
# save("LH2Drepmeandiff129-4.jld","data",repmeandiff)
