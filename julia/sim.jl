using Random
using Printf
using ArgParse
using Distributions

const cell_type = Float64
const cell_update = Tuple{Int, Int, cell_type}

function main(argv)
    nrep = 1
    size = 100
    p = 0.1
    mutsize = 0.1
    specrate = 0.0001
    invrate = 0.2
    timescale = 100 * size / p
    nsteps = 100
    endtime = timescale / nsteps
    outfile = ""
    start_time = time()

    options = ArgParseSettings()
    @add_arg_table options begin
        "--size"
            arg_type = Int
            help = "length of one side of the grid"
            default = 100
        "--specrate"
            help = "speciation rate"
            arg_type = Float64
            default = 0.0001
        "--reps"
            help = "number of repetitions"
            arg_type = Int
            default = 1
        "--mutsize"
            help = "maximum change in mutation event"
            arg_type = Float64
            default = 0.1
        "--outfile"
            help = "output file name"
            arg_type = String
            default = "sim"
    end

    args = ArgParse.parse_args(argv, options)

    if haskey(args, "help")
        println(ArgParse.format_help(options))
        return
    end

    outfile = args["outfile"]
    size = args["size"]
    nrep = args["reps"]
    mutsize = args["mutsize"]
    specrate = args["specrate"]

    timescale = 100 * (size * 1.0 / p)
    endtime = timescale / nsteps

    land_grid_data = Vector{cell_type}(undef, size * size)
    land_grid = [land_grid_data[i*size+1:(i+1)*size] for i in 0:size-1]

    land_mask_data = zeros(size, size)
    land_mask = [land_mask_data[i*size+1:(i+1)*size] for i in 0:size-1]

    println("Inputs:")
    println("\trepetitions = ", nrep)
    println("\tsize = ", size, "x", size)
    println("\tindividuals per patch = ", 1 / p)
    println("\tmutation size = ", mutsize)
    println("\tspeciation rate = ", @sprintf("%.2e", specrate))
    println("\ttimescale = ", timescale)
    println("")

    RanGen = Random.MersenneTwister(1234)

    speciation_distribution = Binomial(size*size, specrate)

    for rep in 1:nrep
        rep_start_time = time()

        fill!(land_grid_data, 1)
        fill!(land_mask_data, 0)

        neighborhood = Vector{cell_type}(undef, 8)
        inv = Vector{cell_type}(undef, 8)
        inv_sum = 0
        inv_index = 0
        i, j, row, col = 0, 0, 0, 0
        grid_min = 1
        grid_max = 1

        for step in 1:timescale
            speciation_event_count = rand(speciation_distribution)
            spec_events = Set{Int}()
            
            while length(spec_events) < speciation_event_count
                index = rand(1:size*size)
                if !(index in spec_events)
                    i = (index - 1) % size + 1
                    j = div(index - 1, size) + 1
                    push!(spec_events, index)
                    down = rand(Bool)
                    ratio = (1 + rand() * mutsize)
                    
                    if down
                        ratio = 1 / ratio
                    end
                    
                    probsuccess = p * ratio / (p * (ratio - 1) + 1)
                    
                    if rand() ≤ probsuccess
                        land_grid[i][j] *= ratio
                        
                        if land_grid[i][j] > grid_max
                            grid_max = land_grid[i][j]
                        end
                        
                        if land_grid[i][j] < grid_min
                            grid_min = land_grid[i][j]
                        end
                        
                        for x in -1:1, y in -1:1
                            row = i + x
                            col = j + y
                            
                            if row > size || row < 1
                                row = (row + size - 1) % size + 1
                            end
                            
                            if col > size || col < 1
                                col = (col + size - 1) % size + 1
                            end
                            
                            local_max = land_grid[row][col]
                            
                            for xx in -1:1, yy in -1:1
                                local_max = max(local_max, land_grid[(row+xx+size-1)%size+1][(col+yy+size-1)%size+1])
                            end
                            
                            land_mask[row][col] = local_max * p / (local_max * p + land_grid[row][col] * (1-p))
                        end
                    end
                end
            end

            for i in 1:size, j in 1:size
                if land_mask[i][j] ≠ 0
                    randval = rand()
                    
                    if randval < land_mask[i][j]
                        inv_sum = 0
                        inv_index = 0
                        
                        for x in -1:1, y in -1:1
                            if (x ≠ 0 || y ≠ 0)
                                row = i + x
                                col = j + y
                                
                                if row > size || row < 1
                                    row = (row + size - 1) % size + 1
                                end
                                
                                if col > size || col < 1
                                    col = (col + size - 1) % size + 1
                                end
                                
                                neighborhood[inv_index+1] = land_grid[row][col]
                                inv[inv_index+1] = p * neighborhood[inv_index+1] / (p * neighborhood[inv_index+1] + land_grid[i][j] * (1-p))
                                inv_sum += inv[inv_index+1]
                                inv_index += 1
                            end
                        end
                        
                        if randval ≤ inv_sum / 8
                            weighted_rand = rand() * inv_sum
                            inv_index = 1
                            
                            while weighted_rand > inv[inv_index]
                                weighted_rand -= inv[inv_index]
                                inv_index += 1
                            end
                            
                            if neighborhood[inv_index] ≠ land_grid[i][j]
                                land_grid[i][j] = neighborhood[inv_index]
                            end
                        end
                    end
                end
            end

            if maximum(land_grid_data) / (1 + mutsize) > typemax(Float64) # Check for overflow
                land_grid_mean = sum(land_grid_data) / (size * size)
                grid_max = 1
                grid_min = 1
                
                for i in 1:size, j in 1:size
                    land_grid[i][j] /= land_grid_mean
                    
                    if land_grid[i][j] > grid_max
                        grid_max = land_grid[i][j]
                    end
                    
                    if land_grid[i][j] < grid_min
                        grid_min = land_grid[i][j]
                    end
                end
                
                println("Normalized at step ", step)
            end
            
            if step % endtime == 0
                println("Progress: ", step / endtime, "%")
            end
        end

        rep_end_time = time()
        rep_time = rep_end_time - rep_start_time
        println("Rep ", rep, " run time = ", rep_time, "s")

        out_start_time = time()
        outfile = args["outfile"] * "_rep$rep.csv"
        open(outfile, "w") do file
            for i in 1:size
                for j in 1:size-1
                    print(file, land_grid[i][j])
                    print(file, ", ")
                end
                println(file, land_grid[i][size])
            end
        end
        out_end_time = time()
        out_time = out_end_time - out_start_time
        println("Output run time = ", out_time, "s")
    end

    end_time = time()
    total_time = end_time - start_time
    println("Total run time = ", total_time, "s")
    println("Wrote results to file ", outfile)
end

main(ARGS)
