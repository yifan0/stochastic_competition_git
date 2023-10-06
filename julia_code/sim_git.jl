# Simulate 2D life history ecosystem with speciation and invasion rules
using DynamicGrids, Crayons, StatsBase, Plots, DataFrames, GLM, Random, JLD
Random.seed!(1234)
# 10 independent simulations of 129X129 system
# no biological reason for using 10
# size should be large enough (compared with speciation rate) such that there are enough species generated.
nrep=1
size = 100
# Each patch is initialized to have a species with 1.0 effective population (fitness)
landrept = fill(1.0,size,size,nrep)
# 1/p individuals in each patch
# p could be other values, e.g., p = 0.5. Changing p affects the mutation & invasion rates
p = .1
# mutation size is 0.1
# no biological reason for 0.1. smaller mutation size may give less fluctuating results in variance measurement
mutsize = 0.1
# speciation rate is 1.0e-4
# We try to push speciation rate to be as small as possible
specrate = 1.0e-4

# timescale of simulation is proportional to size
# According to KPZ results, theoretically the system converges after such time period.
# According to our simulation, the system has already converged in 1/10 of the timescale below.
timescale = floor(Int, 100*(1/p)*size^(1))
# divide the total simulation time into 100 steps. We do renormalization after each step so that the effective population values are not too high.
# no biological reason for 100
nsteps=100
# simulation time of each step
endtime = floor(Int,timescale/nsteps)
println(timescale,endtime)


# Speciation rule (uses DynamicGrids package)
spec_rule = let prob_mut=specrate, mutsize = 0.1
    Cell() do data, cell, I
        randvalue = rand()
        # a mutation occurs with probability prob_mut*probsuccess, so if the randvalue > prob_mut, the mutation cannot occur
        if randvalue < prob_mut
        # Effecive population could either increase or decrease in a speciation event
        upordown=sample([1,-1])
        # The ratio of effective populatoin between new species and old one
        ratio = (1+rand()*mutsize)^upordown
        # Probability that the new species outcompetes the old one:
        # new_eff_pop * p / (new_eff_pop * p + old_eff_pop * (1-p)). After simplification we get:
        probsuccess = p*ratio/(p*(ratio-1)+1)
        # If the new species wins, replace the cell with the new effective population, otherwise keep the original value
        randvalue <= prob_mut*probsuccess ? cell*ratio : cell
        else
            cell
        end
end
end

# Invasion rule (uses DynamicGrids package)
nbr_rule = let prob_inv=1
    Neighbors(Moore(1)) do data, neighborhood, cell, I
        randvalue = rand()
        # an invasion occurs when randvalue < mean(inv)
        # suppose neighbor values are not 2X larger than the current cell, then inv calculated below is < 0.18, so if randvalue > 0.2, the invasion cannot occur.
        if randvalue < 0.2
            # Label neighbors
            nbrvalues = [neighborhood[1],neighborhood[2],neighborhood[3],neighborhood[4],neighborhood[5],neighborhood[6],neighborhood[7],neighborhood[8]]
            # neighbor_eff_pop * p / (neighbor_eff_pop * p + self_eff_pop * (1-p)).
            inv=p*nbrvalues./(p*nbrvalues+fill(cell*(1-p),8)) 
            # If invasion happens, replace with the neighbor's effective population
            if randvalue <= mean(inv)
                # Probability of invasion labeled with neighbor index
                invindex = wsample(collect(1:8),inv)
                neighborhood[invindex]
            else
                cell
            end
        else
            cell
        end
end
end



function simulation(nrep,step)
    # Do simulations for each copy at specific step
    # Multi-threading
    Threads.@threads for n in 1:nrep
        # Initiate as the nth copy
        init = landrept[:,:,n]
        # Save the simulation result at the endtime
        output = ResultOutput(init; tspan=1:endtime)
        # Simulates with speciation and invasion rules, uses ThreadedCPU, set the boundary to be wrapped (sim! function is defined in DynamicGrids package)
        sim!(output, spec_rule,nbr_rule; proc = ThreadedCPU(), boundary=Wrap())
        init=output[1]
        # renormalize the simulated result so that the fittest species will not have a very high effective population so that it always outcompetes other species.
        init=init/mean(init)
        # replace the nth copy with the simulated result
        landrept[:,:,n] = init
    end
    println(step)
end

# Wrap into a single function for performance
function func(nsteps, nrep)
    # Do simulations for each steps
    for step in 1:nsteps
        simulation(nrep,step)
    end
end


@time begin
    func(nsteps, nrep)
end

# Save the 10 independent copies of 129*129 system
save("LH2D129-4.jld", "data",landrept)


