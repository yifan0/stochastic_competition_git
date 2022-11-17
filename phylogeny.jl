using DynamicGrids, Crayons, StatsBase, Plots, DataFrames, GLM, Random, JLD
Random.seed!(1234)
# simulate 10 times
# nrep=10
# divide the total simulation time into 100 steps
nsteps=100
# 129*129 system
size = 129
p = .1
# mutation size is 0.1
mutsize = 0.1
# speciation rate is 1.0e-4
specrate = 1.0e-3
# timescale is proportional to size
timescale = floor(Int, 100*(1/p)*size^(1))
#
endtime = floor(Int,timescale/nsteps)
println(timescale,endtime)

landscape =Array{Float64,2}(undef,size,size)
landrept = fill(1.0,size,size,1)


phylodata = Array{Any,2}(missing,0,3)
tempori=zeros(0)
tempmut=zeros(0)

# speciation rule
spec_rule = let prob_mut=specrate, mutsize = 0.1
    Cell() do data, cell, I
        upordown=sample([1,-1])
        ratio = (1+rand()*mutsize)^upordown
        probsuccess = p*ratio/(p*(ratio-1)+1)
        if rand() <= prob_mut*probsuccess
            push!(tempori,cell)
            push!(tempmut,cell*ratio)
            cell*ratio
        else
            cell
        end
        # rand() <= prob_mut*probsuccess ? cell*ratio : cell
    end
end

# invasion rule
nbr_rule = let prob_inv=1
    Neighbors(Moore(1)) do data, neighborhood, cell, I
        nbrvalues = [neighborhood[1],neighborhood[2],neighborhood[3],neighborhood[4],neighborhood[5],neighborhood[6],neighborhood[7],neighborhood[8]]
        inv=p*nbrvalues./(p*nbrvalues+fill(cell*(1-p),8)) 
        invindex = wsample(collect(1:8),inv)
        rand() <= mean(inv) ? neighborhood[invindex] : cell
end
end

function phylo(output,step)
    global avg = mean(output[endtime])
    survived = unique(output[endtime])
    for sur in survived
        # "speciation" by renormalization when loading to next step
        global phylodata = vcat(phylodata,[sur,sur/avg,step*endtime]')
        # record true speciations
        # while the current value does not appear in the first time unit, trace its ancestor
        while !(sur in output[1])
            for t in 1:endtime
                if sur in output[t]
                    index = findfirst(x->x==sur,tempmut)
                    ori = tempori[index]
                    # record ancestor, descendant, time
                    global phylodata = vcat(phylodata,[ori,sur,(step-1)*endtime+t]')
                    # let the ancestor be the descendant and repeat
                    sur = ori
                    break
                end
            end
        end
    end
end


function simulation(step)
        # Threads.@spawn begin
        init = landrept[:,:,1]
        init=init/mean(init)
        output = ArrayOutput(init; tspan=1:endtime)
        sim!(output, spec_rule,nbr_rule; boundary=Wrap())
        phylo(output,step)
        # global out = output
        init=output[endtime]
        landrept[:,:,1] = init
        # end
    println(step)
end

function func(nsteps)
    for step in 1:nsteps
        simulation(step)
        # reset tempmut,tempori
        global tempmut=zeros(0)
        global tempori=zeros(0)
    end
end



@time begin
    func(nsteps)
end

# save("phylo129-3.jld","data",phylodata)
# phylodata = load("phylo129-3.jld")["data"]

# now we have a dataframe of ancestor; descendant; time of speciation 
df = DataFrame(original=phylodata[:,1],mutation=phylodata[:,2],time=phylodata[:,3])
df = sort(df,:time)
df = unique!(select(df,[:original,:mutation,:time]))

# create a history tree for each survived species
function genTrees(species)
    nspecies = length(species)
    trees = Array{Any,1}(missing,nspecies)
    for n in 1:nspecies
        # current species
        current = species[n]
        trees[n] = DataFrame(original=Float64[],mutation=Float64[],time=Float64[])
        while (current != 1.0)
            # trace back
            row = df[df.mutation.==current,:]
            trees[n] = vcat(trees[n],row)
            # replace
            current = row[1,1]
        end
    end
    return trees
end

# find common ancestor between species m and species n so that know when and how to merge
function ancestor(m,n,ancestorT,trees)
    for i in 1:length(trees[m][:,1])
        for j in 1:length(trees[n][:,1])
            # for the case that the species itself is the common ancestor
            if trees[m][i,2] == trees[n][j,2]
                anc = trees[m][i,2]
                # first occurence of common ancestor in each tree
                indexm = findfirst(x->x==anc,trees[m][:,2])
                indexn = findfirst(x->x==anc,trees[n][:,2])
                if indexm == 1
                    tm = timescale
                else
                    tm = trees[m][indexm-1,3]
                end
                if indexn == 1
                    tn = timescale
                else
                    tn = trees[n][indexn-1,3]
                end
                # time of occurence
                return ancestorT[m,n,1],ancestorT[m,n,2] = anc, min(tm,tn)
            end
        end
    end
end


# merge closest species
@time begin
# all survived species
species = df[df.time.==timescale,:][:,1]
# count the number of patches occupied by each species
proportions = [count(x->x==i,landrept) for i in species]
# attach number of patches occupied to corresponding species
code = [string(i) for i in proportions]
# dictionary mapping  species to their names
dict = Dict(zip(species,code))
# dictionary mapping species to speciation time
dictTime = Dict(zip(species,[timescale for i in species]))
# the newick format phylogeny tree we need
global lastcode = Array{String,1}(undef,0)
# merge when there are more than one species
while (length(species) > 1)
println(length(species))
ancestorT = Array{Any,3}(missing,length(species),length(species),2)
trees = genTrees(species)
# double loop iterate each pair
Threads.@threads for m in 1:length(species)-1
    Threads.@threads for n in m+1:length(species)
        # find common ancestor for each pair
        ancestor(m,n,ancestorT,trees)
    end
end
# find the latest merge
maxT,maxCoord = findmax(skipmissing(ancestorT[:,:,2]))
maxAnc = ancestorT[maxCoord[1],maxCoord[2],1]
toMergeIndex = findall(t->t==maxT,skipmissing(ancestorT[:,:,2]))
toMergeIndexUnique = sort(unique([toMergeIndex[i][j] for i in eachindex(toMergeIndex), j in 1:2]))
# time difference is the edge length in phylogeny tree 
timeDiff = [get(dictTime,species[i],0) for i in toMergeIndexUnique].-maxT
# name the common ancestor
toMerge = [get(dict,species[toMergeIndexUnique[i]],0)*":"*string(timeDiff[i]) for i in eachindex(timeDiff)]
newcode = "("*join(string.(toMerge),",")*")"
# add the common ancestor
push!(dict,maxAnc=>newcode)
push!(dictTime,maxAnc=>maxT)
# remove the merged species
deleteat!(species,toMergeIndexUnique)
push!(species,maxAnc)
# end when all species are traced back to one common ancestor 
if length(species) == 1
    lastcode = newcode*";"
end
end

end

open("phylo129-3.txt","w") do file
    write(file,lastcode)
end
