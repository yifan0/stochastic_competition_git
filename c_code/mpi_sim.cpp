#include <cstdio>
#include <random> // rand()
#include <string.h> // memcpy()
#include <queue> // queue
#include <omp.h>
#include <time.h>
#include <chrono>
#include <set>
#include <map>
#include <numeric>
#include <algorithm>
#include "sfmt.h"
#include "sfmt.cpp"
#include "userintf.cpp"
#include <mpi.h>
#include <iostream> // std::cout
#include <cxxopts.hpp> // to handle cmdline args
#include "SimGrid.cpp"
using namespace std;

#define println(...) { if(rank == 0) { printf(__VA_ARGS__); printf("\n"); } }
#define debug(...) { if(true) { printf("Rank %d: ", rank); printf(__VA_ARGS__); printf("\n"); } }
#define print(...) { printf(__VA_ARGS__); }

#define   NDIM         2
#define   GHOSTS       2

typedef double cell_type;
typedef std::tuple<int, int, cell_type> cell_update;

int main(int argc, char* argv[]){
    int nrep = 10;
    int size = 100;
    double p = 0.1;
    double mutsize = 0.1;
    double specrate = 0.0001; // 1.0e-4
    double invrate = 0.2;
    int timescale = 100*size/p;
    int nsteps = 100;
    int endtime = timescale/nsteps;
    std::string outfile;
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time, rep_start_time, rep_end_time, out_start_time, out_end_time;
    int rank, nprocs;
    FILE *fp;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cxxopts::Options options("mpirun -n N ./ga_sim", "Stocastic competition simulation with MPI parallelism");

    options.add_options()
        ("s,size", "grid size", cxxopts::value<int>()->default_value("500"))
        ("r,reps", "number of repetitions", cxxopts::value<int>()->default_value("1"))
        ("o,outfile", "output file prefix", cxxopts::value<std::string>()->default_value("out"))
        ("h,help", "Print usage")
        ;

    auto result = options.parse(argc, argv);

    if(result.count("help")) {
        if(rank == 0)
            std::cout << options.help() << std::endl;
        MPI_Finalize();
        exit(0);
    }
    size = result["size"].as<int>();
    nrep = result["reps"].as<int>();

    timescale = 100*(size*1.0/p);
    endtime = timescale/nsteps;

    start_time = std::chrono::system_clock::now();

    SimGrid<cell_type> grid(size, GHOSTS);
    SimGrid<bool> mask(size, 0);

    println("Inputs:")
        println("\trepetitions = %d", nrep)
        println("\tsize = %d%s%d", size, "x", size)
        println("\tindividuals per patch = %f", 1/p)
        println("\tmutation size = %f", mutsize)
        println("\tspeciation rate = %f", specrate)
        println("\tinvasion rate = %f", invrate)
        println("\ttimescale = %d", timescale)
        println("\tnsteps = %d", nsteps)
        println("\tend time = %d", endtime);
    println("");
    fflush(stdout);

    // Random number generation
    CRandomSFMT1 RanGen(time(0)+rank*10); // Agner Combined generator
    // distribution for speciation events
    std::default_random_engine generator;
    std::binomial_distribution<int> speciation_distribution(size*size, specrate);

    for(int rep = 0; rep < nrep; rep++) {
        rep_start_time = std::chrono::system_clock::now();

        // TODO: get rows and columns from the grids
        debug("rows = %d, cols = %d", grid.getRows(), grid.getCols());
        for(size_t i = 0; i < grid.getRows(); i++) {
            for(size_t j = 0; j < grid.getCols(); j++) {
                grid[i][j] = 1;
                mask[i][j] = 0;
            }
        }
        debug("filled in default values");

        // invasion rule variables
        cell_type neighborhood[8];
        cell_type inv[8];
        cell_type inv_sum = 0;
        int inv_index = 0;

        for(int step = 0; step < timescale; step++) {

            // speciation rule
            int speciation_event_count = speciation_distribution(generator);
            std::set<int> spec_events;
            while(spec_events.size() < speciation_event_count) {
                int index = rand() % ((grid.getRows())*(grid.getCols()));
                if(!spec_events.count(index)) {
                    int i = index % (grid.getRows());
                    int j = index / (grid.getCols());
                    spec_events.insert(index);
                    int down = rand() % 2;
                    float ratio = (1+RanGen.Random()*mutsize);
                    if(down) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(RanGen.Random() <= probsuccess) {
                        grid[i][j] *= ratio;
                        mask[i][j] = true;
                        for(int x = -1; x <= 1; x++) {
                            for(int y = -1; y <= 1; y++) {
                                if((x != 0 || y != 0)) {
                                    mask[i+x][j+y] = true;
                                }
                            }
                        }
                    }
                }
            }

            if(step%10 == 0) {
                grid.sync();
                mask.sync();
            }

            // invasion rule
            std::vector<cell_update> updates;
            for(int i = 0; i < grid.getRows(); i++) {
                for(int j = 0; j < grid.getCols(); j++) {
                    if(mask[i][j] || i == 0 || i == mask.getRows() || j == 0 || j == mask.getCols()) {
                        double randval = RanGen.Random(); //random_float();
                        if(randval < invrate) {
                            inv_sum = 0;
                            inv_index = 0;
                            for(int x = -1; x <= 1; x++) {
                                for(int y = -1; y <= 1; y++) {
                                    if((x != 0 || y != 0)) {

                                        neighborhood[inv_index] = grid[i+x][j+y];
                                        inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+grid[i][j]*(1-p));
                                        inv_sum += inv[inv_index];
                                        inv_index++;
                                    }
                                }
                            }

                            if(randval <= inv_sum/8) {
                                // Get random element with weighted probabilities
                                double weighted_rand = RanGen.Random()*inv_sum;
                                inv_index = 0;
                                while(weighted_rand > inv[inv_index]) {
                                    weighted_rand -= inv[inv_index];
                                    inv_index++;
                                }
                                if(neighborhood[inv_index] != grid[i][j]) {
                                    updates.push_back({i, j, neighborhood[inv_index]});
                                }
                            }
                        }
                    }
                }
            }

            for(const auto& [i, j, val] : updates) {
                grid[i][j] = val;
                bool unmask = true;
                for(int x = -1; x <= 1; x++) {
                    for(int y = -1; y <= 1; y++) {
                        if((x != 0 || y != 0)) {
                            bool unmask = true;
                            for(int xx = -1; xx <= 1; xx++) {
                                for(int yy = -1; yy <= 1; yy++) {
                                    if((xx != 0 || yy != 0)) {
                                        if(grid[i+x+xx][j+y+yy] != val) {
                                            unmask = false;
                                            break;
                                        }
                                    }
                                }
                                if(!unmask) break;
                            }
                            mask[i+x][j+y] = !unmask;
                        }
                    }
                }
            }

            // renormalize every nstep steps
            if(step%endtime == endtime-1) {
                cell_type land_grid_mean = 0;
                for(int i = 0; i < grid.getRows(); i++) {
                    for(int j = 0; j < grid.getCols(); j++) {
                        land_grid_mean += grid[i][j];
                    }
                }
                land_grid_mean /= (size*size);
                MPI_Allreduce(MPI_IN_PLACE, &land_grid_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for(int i = 0; i < grid.getRows(); i++) {
                    for(int j = 0; j < grid.getCols(); j++) {
                        grid[i][j] /= land_grid_mean;
                    }
                }
            }
        }

        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);

        out_start_time = std::chrono::system_clock::now();
        /*
        outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + ".log";
        fp = fopen(outfile.c_str(), "w");
        GA_Print_file(fp, ga_land_grid);
        fclose(fp);
        */
        out_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> out_time = out_end_time - out_start_time;
        println("Output run time = %fs", out_time.count());

    }

    end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;
    println("Total run time = %fs", total_time.count());
    println("Wrote results to file %s", outfile.c_str());
    fflush(stdout);

    MPI_Finalize();
}

