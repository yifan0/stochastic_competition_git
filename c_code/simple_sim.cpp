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
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

#define USAGE \
"sim gridsize nreps outfile"

#define HELPTEXT \
"Simulate plant spread on grid\n\
    gridsize : number of cells along an edge of the grid\n\
    nreps : number of repetitions of the simulation\n\
    outfile : output file name\n"

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

    bool badargs = false;
    if(argc < 1 || argc > 4) badargs = true;
    if(argc > 1) {
        size = stoi(argv[1]);
        if(size < 1) {
            badargs = true;
        }
    }
    if(argc > 2) {
        nrep = stoi(argv[2]);
        if(nrep < 1) {
            badargs = true;
        }
    }
    if(argc > 3) {
        outfile = argv[3];
    }
    timescale = 100*(size*1.0/p);
    endtime = timescale/nsteps;

    if(badargs) {
        fprintf(stderr, ">E Usage: %s\n", USAGE);
        exit(1);
    }

    start_time = std::chrono::system_clock::now();

    // grid for average across reps
    std::vector<cell_type> mean_grid_data(size*size, 0); // fill grid with 0s
    std::vector<cell_type*> mean_grid_arrays;
    for(int i = 0; i != size*size; i += size) {
        mean_grid_arrays.push_back(mean_grid_data.data() + i);
    }
    cell_type** mean_grid = mean_grid_arrays.data();

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
    CRandomSFMT1 RanGen((int)time(NULL)); // Agner Combined generator
    // distribution for speciation events
    std::default_random_engine generator;
    std::binomial_distribution<int> speciation_distribution(size*size, specrate);

    for(int rep = 0; rep < nrep; rep++) {
        rep_start_time = std::chrono::system_clock::now();
        // grid for land data
        // strategy from: https://stackoverflow.com/questions/46354262/how-to-dynamically-allocate-a-contiguous-2d-array-in-c
        std::vector<cell_type> land_grid_data(size*size);
        std::vector<cell_type*> land_grid_arrays;
        for(int i = 0; i != size*size; i += size) {
            land_grid_arrays.push_back(land_grid_data.data() + i);
        }
        cell_type** land_grid = land_grid_arrays.data();

        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid[i][j] = 1;
            }
        }

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
                int index = rand() % (size*size);
                if(!spec_events.count(index)) {
                    int i = index % size;
                    int j = index / size;
                    spec_events.insert(index);
                    int down = rand() % 2;
                    float ratio = (1+RanGen.Random()*mutsize);
                    if(down) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(RanGen.Random() <= probsuccess) {
                        land_grid[i][j] *= ratio;
                    }
                }
            }

            // invasion rule
            std::vector<cell_update> updates;
            for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++) {
                    double randval = RanGen.Random(); //random_float();
                    if(randval < invrate) {
                        inv_sum = 0;
                        inv_index = 0;
                        for(int x = -1; x < 1; x++) {
                            for(int y = -1; y <= 1; y++) {
                                if((x != 0 || y != 0) && i+x >= 0 && i+x < size && j+y >= 0 && j+y < size) {

                                    neighborhood[inv_index] = land_grid[i+x][j+y];
                                    inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+land_grid[i][j]*(1-p));
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
                            if(neighborhood[inv_index] != land_grid[i][j]) {
                                //land_grid[i][j] = neighborhood[inv_index];
                                updates.push_back({i, j, neighborhood[inv_index]});
                            }
                        }
                    }
                }
            }

            for(const auto& [i, j, val] : updates) {
                land_grid[i][j] = val;
            }

            // renormalize every nstep steps
            if(step%nsteps == 0) {
                float land_grid_mean = 0;
                float tmp_sum = 0;
                for(int i = 0; i < size; i++) {
                    for(int j = 0; j < size; j++) {
                        tmp_sum += land_grid[i][j];
                    }
                    land_grid_mean += tmp_sum/(size*size);
                    tmp_sum = 0;
                }
                for(int i = 0; i < size; i++) {
                    for(int j = 0; j < size; j++) {
                        land_grid[i][j] /= land_grid_mean; // normalize grid
                        mean_grid[i][j] += land_grid[i][j]/nrep; // add contribution to average across reps
                    }
                }
            }

        }

        // normalize the simulated result so that the fittest species will not have a very high effective population so that it always outcompetes other species
        float land_grid_mean = 0;
        float tmp_sum = 0;
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                tmp_sum += land_grid[i][j];
            }
            land_grid_mean += tmp_sum/(size*size);
            tmp_sum = 0;
        }
        println("Grid mean = %f for rep %d", land_grid_mean, rep);
        fflush(stdout);
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid[i][j] /= land_grid_mean; // normalize grid
                mean_grid[i][j] += land_grid[i][j]/nrep; // add contribution to average across reps
            }
        }
        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);
    }

    out_start_time = std::chrono::system_clock::now();
    FILE *fp;
    fp = fopen(outfile.c_str(), "w");
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size-1; j++) {
            fprintf(fp, "%f", mean_grid[i][j]);
            fprintf(fp, ", ");
        }
        fprintf(fp, "%f\n", mean_grid[i][size-1]);
    }
    fclose(fp);
    out_end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> out_time = out_end_time - out_start_time;
    println("Output run time = %fs", out_time.count());
    fflush(stdout);
    end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;
    println("Total run time = %fs", total_time.count());
    println("Wrote reuslts to file %s", outfile.c_str());
    fflush(stdout);
}

