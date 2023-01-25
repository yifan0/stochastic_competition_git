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

enum DIRECTION { NW=0, N=1, NE=2, W=10, HERE=11, E=12, SW=20, S=21, SE=22 };

typedef double cell_type;

int ns_dir(DIRECTION dir) {
    if (dir % 3 == 0) return -1;
    if (dir % 3 == 1) return 0;
    if (dir % 3 == 2) return 1;
}

int ew_dir(DIRECTION dir) {
    if (dir % 3 == 0) return -1;
    if (dir % 3 == 10) return 0;
    if (dir % 3 == 20) return 1;
}

struct cellUpdateRecord {
    int x;
    int y;
    cell_type val;
};

cellUpdateRecord newCellUpdateRecord(int x, int y, cell_type val) {
    cellUpdateRecord rec;
    rec.x = x;
    rec.y = y;
    rec.val = val;
    return rec;
}

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
    if(argc < 3 || argc > 4) badargs = true;
    if (argc < 2 || sscanf(argv[1],"%ld",&size) != 1 || size < 1)
        badargs = true;
    if (argc < 3 || sscanf(argv[2],"%ld",&nrep) != 1 || nrep < 1)
        badargs = true;
    if(argc > 3)
       outfile = argv[3];
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
    // distribution for speciation events
    std::default_random_engine generator;
    std::binomial_distribution<int> speciation_distribution(size*size, specrate);
    std::binomial_distribution<int> invasion_distribution(size*size, 0.2);
    std::uniform_real_distribution<double> real_distribution(0.0, 1.0);
    std::uniform_real_distribution<double> speciation_real_distribution(0.0, 0.2);

    // Random generator for shuffle
    std::random_device rd;
    std::mt19937 g(rd());
    //std::minstd_rand0 minrand;
    //boost::random::mt19937 dis;
    //boost::random::uniform_real_distribution<double> gen(0.0, 1.0);
    //std::random_device rd;
    //std::mt19937 gen(rd());
    //std::uniform_real_distribution<> dis(0.0, 1.0);
    
    // List of all indices to select from
    std::vector<unsigned> all_cell_indices(size*size);
    std::iota(all_cell_indices.begin(), all_cell_indices.end(), 0);
    
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
                    float ratio = (1+real_distribution(generator)*mutsize);
                    if(down) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(speciation_real_distribution(generator) < specrate*probsuccess) { // random value must be < 0.2 to have made it here
                        land_grid[i][j] *= ratio;
                    }
                }
            }
           
            // invasion rule
            //std::map<int, cell_type> invasion_events;
            for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++) {
                    //double randval = 0.4;
                    double randval = real_distribution(generator);
                    if(false && randval < invrate) {
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
                            std::uniform_real_distribution<double> weighted_dist(0, inv_sum);
                            double weighted_rand = weighted_dist(generator);
                            inv_index = 0;
                            while(weighted_rand > inv[inv_index]) {
                                weighted_rand -= inv[inv_index];
                                inv_index++;
                            }
                            land_grid[i][j] = neighborhood[inv_index];
                            // TODO: check for references to already changed cells -- instead of saving the modifications, try saving the originals and checking for conflicts
                            //invasion_events[i*size + j] = neighborhood[inv_index];
                        }
                    }
                }
            }

            // For every invasion event, update the grid
            //println("Number of invasion events on step %d = %d", step, invasion_events.size());
            //fflush(stdout);
            /*
            for(auto const& event : invasion_events) {
                int i = event.first % size;
                int j = event.first / size;
                land_grid[i][j] = event.second;
            }
            */

            /*
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
            */

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

