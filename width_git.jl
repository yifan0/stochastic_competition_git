# Measure diversity & width at sample regions
using DynamicGrids, Crayons, StatsBase, Plots, DataFrames, GLM, Random, JLD
Random.seed!(1234)
# 10 independent copies of 129X129 system
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
    # Measure diversity/width of the sample region, store as the (j1,j2,i) element of diversity/width matrix
    diversity[j1,j2,i]=length(unique(@view(landscape[srange1,srange2])))-1
    width[j1,j2,i] = var(@view(landscape[srange1,srange2]))*(i^2+2*i)/(i+1)^2
end
# Average over all size i samples
smeandiv[i,n]=mean(skipmissing(@view(diversity[:,:,i])))
smeanwidth[i,n]=mean(skipmissing(@view(width[:,:,i])))
end
# Average over all copies of simulations
repmeandiv[i]= mean(@view(smeandiv[i,:]))
repmeanwidth[i]= mean(@view(smeanwidth[i,:]))
end

# Save data
# save("LH2Drepmeandiv129-4.jld","data",repmeandiv)
# save("LH2Drepmeanwidth129-4.jld","data",repmeanwidth)