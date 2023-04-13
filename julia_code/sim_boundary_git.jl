using Agents, StatsBase, Random
size = 128
nthreads = 8
nrep = 10
# slabsize = size*Int(size/nthreads)
p = 0.1
mutsize = 0.1
specrate = 1.0e-4
timescale = floor(Int,100*(1/p)*size)
nsteps = 100
endtime = floor(Int,timescale/nsteps)
@agent SpeciesAgent GridAgent{2} begin
    value::Float64 # effective population
    status::Bool # on "boundary" or not
end

# @everywhere mutable struct SpeciesAgent <: AbstractAgent
#     id::Int             
#     pos::NTuple{2, Int} 
#     value::Float64
# end

function initialize(; numagents = size*size, griddims = (size,size), seed = 125)
    space = GridSpaceSingle(griddims; periodic = false)
    rng = Random.MersenneTwister(seed)
    model = ABM(SpeciesAgent, space; rng)
    for n in 1:numagents
        agent = SpeciesAgent(n, (mod1(n,size),ceil(Int,n/size)), 1.0, false)
        add_agent_pos!(agent, model)
    end
    return model
end

function agent_step!(agent,model)
    speciation!(agent,model)
    invasion!(agent,model)
end

function speciation!(agent,model)
    randvalue = rand()
    if randvalue < specrate
    upordown=sample([1,-1])
    ratio = (1+rand()*mutsize)^upordown
    probsuccess = p*ratio/(p*(ratio-1)+1)
    if randvalue < specrate*probsuccess
        agent.value *= ratio
        # update neighbor boundary status
        for nbr in nearby_agents(agent,model)
            nbr.status = true
        end
    end
    end
end

# check boundary status
function check_boundary!(agent,model)
    nbrvalues = [nbr.value for nbr in nearby_agents(agent,model)]
    for nbrvalue in nbrvalues
        if agent.value != nbrvalue
            agent.status = true
            return
        end
    end
    agent.status = false
end

# Note that when periodic = false, the number of neighbors is not always 8
# assuming neighbor values are not 2X smaller
function invasion!(agent,model)
    # do nothing if not on boundary
    agent.status == false && return
    randvalue = rand()
    if randvalue < 0.2
    nbrvalues = [nbr.value for nbr in nearby_agents(agent,model)]
    nNeighbor = length(nbrvalues)
    inv=p*nbrvalues./(p*nbrvalues+fill(agent.value*(1-p),nNeighbor)) 
        if randvalue <= mean(inv)
            invindex = wsample(collect(1:nNeighbor),inv)
            agent.value = nbrvalues[invindex]
            # update neighbor status
            check_boundary!(agent,model)
            for nbr in nearby_agents(agent,model)
                check_boundary!(nbr,model)
            end
        end
    end
end

# adata = [:pos,:value]
landrept = Array{Float64,3}(undef, size, size, nrep)
models = [initialize() for i in 1:nrep]
@time begin
for step in 1:nsteps
    Threads.@threads for n in 1:nrep
        run!(models[n], agent_step!, endtime)
        landrept[:,:,n] = [models[n][m].value for m in 1:size*size]
        local m = mean(landrept[:,:,n])
        for i in 1:size*size
            models[n][i].value /= m
        end
    end
    println(step)
end
end


# 128X128 10rep 1872s