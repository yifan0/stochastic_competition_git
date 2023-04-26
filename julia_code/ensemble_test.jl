using Distributed
addprocs()
@everywhere using Agents, StatsBase, Random, JLD
@everywhere size = 33
@everywhere p = 0.1
@everywhere mutsize = 0.1
@everywhere specrate = 1.0e-3
@everywhere timescale = floor(Int,100*(1/p)*size)
@everywhere nrep = 8
@everywhere nsteps = 100
@everywhere endtime = floor(Int,timescale/nsteps)

@everywhere @agent SpeciesAgent GridAgent{2} begin
    value::Float64
end

# @everywhere mutable struct SpeciesAgent <: AbstractAgent
#     id::Int             
#     pos::NTuple{2, Int} 
#     value::Float64
# end

@everywhere function initialize(; numagents = size*size, griddims = (size,size), seed = 125)
    space = GridSpaceSingle(griddims, periodic=true)
    rng = Random.MersenneTwister(seed)
    model = ABM(SpeciesAgent, space; rng)
    for n in 1:numagents
        agent = SpeciesAgent(n, (1, 1), 1.0)
        add_agent_single!(agent, model)
    end
    return model
end

#step function
@everywhere function agent_step!(agent,model)
    speciation!(agent)
    invasion!(agent,model)
end

# speciation
@everywhere function speciation!(agent)
    upordown=sample([1,-1])
    ratio = (1+rand()*mutsize)^upordown
    probsuccess = p*ratio/(p*(ratio-1)+1)
    if rand() <= specrate*probsuccess
        agent.value *= ratio
    end
end

# invasion
@everywhere function invasion!(agent,model)
    nbrvalues = [nbr.value for nbr in nearby_agents(agent,model)]
    inv=p*nbrvalues./(p*nbrvalues+fill(agent.value*(1-p),8)) 
    invindex = wsample(collect(1:8),inv)
    if rand() <= mean(inv)
        agent.value = nbrvalues[invindex]
    end
end

@everywhere data = []
@everywhere adata = [:pos,:value]
@everywhere models = [initialize(seed = x) for x in rand(UInt8, nrep)];
@time begin
for step in 1:nsteps    
    global data,sth,models = ensemblerun!(models, agent_step!, dummystep, endtime; adata,parallel = true)
    # ensemblerun!(models, agent_step!, dummystep, endtime; parallel = true)
    for n in 1:nrep
        landscape = [models[n][m].value for m in 1:size*size]
        global m = mean(landscape)
        for i in 1:size*size
            models[n][i].value /= m
        end
    end
    println(step)
end
end

# extract data
# landrept = Array{Float64,3}(undef,size,size,nrep)
# landscape2 = Array{Float64,2}(undef,size,size)
# for n in 1:nrep
#     for id in 1:size*size
#     agent = models[n][id]
#     landscape2[agent.pos[1],agent.pos[2]] = agent.value
#     end
#     landrept[:,:,n] = landscape2
# end

# save("test.jld","data",landrept)