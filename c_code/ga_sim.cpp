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
#include "macdecls.h"
#include <mpi.h>
#include <iostream> // std::cout
#include <cxxopts.hpp> // to handle cmdline args
#include "ga++.h"
using namespace std;

#define println(...) { if(me == 0) { printf(__VA_ARGS__); printf("\n"); } }
//#define debug(...) { if(true) { printf("Rank %d: ", me); printf(__VA_ARGS__); printf("\n"); GA_Sync(); } }
#define debug(...) {}
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
    int me, nprocs;
    int dims[NDIM];
    int lo[NDIM], hi[NDIM];
    int mask_lo[NDIM], mask_hi[NDIM];
    int grid_ld[NDIM], mask_ld[NDIM];
    int ghost_mask_ld[NDIM-1];
    int ghost_grid_ld[NDIM-1];
    int ghost_dims[NDIM], ghost_mask_dims[NDIM];
    int ghost_width[NDIM];
    double* land_grid_ptr, *ghost_grid_ptr;
    int* land_mask_ptr, *ghost_mask_ptr;
    int heap=3000000, stack=3000000; // TODO: find out if these are good values

    MPI_Init(&argc, &argv);
    GA_Initialize();
    if(! MA_init(C_DBL, stack, heap) ) GA_Error("MA_init failed",stack+heap);

    me = GA_Nodeid();
    nprocs = GA_Nnodes();

    cxxopts::Options options("mpirun -n N ./ga_sim", "Stocastic competition simulation with GlobalArray parallelism");

    options.add_options()
        ("s,size", "grid size", cxxopts::value<int>()->default_value("500"))
        ("r,reps", "number of repetitions", cxxopts::value<int>()->default_value("1"))
        ("o,outfile", "output file prefix", cxxopts::value<std::string>()->default_value("out"))
        ("h,help", "Print usage")
        ;

    auto result = options.parse(argc, argv);

    if(result.count("help")) {
        if(me == 0)
            std::cout << options.help() << std::endl;
        MPI_Finalize();
        exit(0);
    }
    outfile = result["outfile"].as<std::string>();
    size = result["size"].as<int>();
    nrep = result["reps"].as<int>();

    timescale = 100*(size*1.0/p);
    endtime = timescale/nsteps;

    start_time = std::chrono::system_clock::now();

    // grid for average across reps
    GA_Mask_sync(1, 1); // turns on/off sync when updating ghosts. TODO: turn off for both before and after and check performance
    for(int i = 0; i < NDIM; i++) {
        dims[i] = size;
        ghost_width[i] = GHOSTS;
    }
    int ga_land_grid = NGA_Create_ghosts(C_DBL, NDIM, dims, ghost_width, "land grid", NULL);
    int ga_land_mask = NGA_Create_ghosts(C_INT, NDIM, dims, ghost_width, "land mask", NULL);

    NGA_Distribution(ga_land_grid, me, lo, hi);
    NGA_Distribution(ga_land_mask, me, mask_lo, mask_hi);
    for(int i = 0; i < NDIM; i++) {
        if(lo[i] != mask_lo[i] || hi[i] != mask_hi[i]) {
            GA_Error("Land grid and mask must have same distribution", 1);
        }
    }

    GA_Sync();
    NGA_Update_ghosts(ga_land_grid);
    NGA_Update_ghosts(ga_land_mask);

    NGA_Access(ga_land_grid, lo, hi, &land_grid_ptr, grid_ld);
    NGA_Access(ga_land_mask, mask_lo, mask_hi, &land_mask_ptr, mask_ld);
    NGA_Access_ghosts(ga_land_grid, ghost_dims, &ghost_grid_ptr, ghost_grid_ld);
    NGA_Access_ghosts(ga_land_mask, ghost_mask_dims, &ghost_mask_ptr, ghost_mask_ld);

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

        int zero = 0;
        double one = 1;
        debug("starting Fill");
        GA_Fill(ga_land_grid, &one);
        GA_Fill(ga_land_mask, &zero);

        // invasion rule variables
        cell_type neighborhood[8];
        cell_type inv[8];
        cell_type inv_sum = 0;
        int inv_index = 0;

        for(int step = 0; step < timescale; step++) {

            debug("Starting step %d", step);

            // speciation rule
            int speciation_event_count = speciation_distribution(generator);
            std::set<int> spec_events;
            while(spec_events.size() < speciation_event_count) {
                int index = rand() % ((hi[0]-lo[0])*(hi[1]-lo[1]));
                if(!spec_events.count(index)) {
                    int i = index % (hi[0]-lo[0]);
                    int j = index / (hi[1]-lo[1]);
                    spec_events.insert(index);
                    int down = rand() % 2;
                    float ratio = (1+RanGen.Random()*mutsize);
                    if(down) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(RanGen.Random() <= probsuccess) {
                        land_grid_ptr[(i+GHOSTS)*grid_ld[0]+j+GHOSTS] *= ratio;
                        land_mask_ptr[(i+GHOSTS)*mask_ld[0]+j+GHOSTS] = true;
                        for(int x = -1; x < 1; x++) {
                            for(int y = -1; y <= 1; y++) {
                                if((x != 0 || y != 0)) {
                                    land_mask_ptr[ (i+x+GHOSTS)*mask_ld[0] + j+y + GHOSTS] = true;
                                }
                            }
                        }
                    }
                }
            }

            debug("After speciation %d", step);

            GA_Update_ghosts(ga_land_grid);
            GA_Update_ghosts(ga_land_mask);

            debug("After update ghosts %d", step);

            // invasion rule
            std::vector<cell_update> updates;
            for(int i = 0; i < hi[0]-lo[0]; i++) {
                for(int j = 0; j < hi[1]-lo[1]; j++) {
                    if(land_mask_ptr[i*mask_ld[0]+j]) {
                        double randval = RanGen.Random(); //random_float();
                        if(randval < invrate) {
                            inv_sum = 0;
                            inv_index = 0;
                            for(int x = -1; x < 1; x++) {
                                for(int y = -1; y <= 1; y++) {
                                    if((x != 0 || y != 0)) {

                                        neighborhood[inv_index] = land_grid_ptr[ (i+x)*grid_ld[0] + j+y ];
                                        inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+land_grid_ptr[i*grid_ld[0]+j]*(1-p));
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
                                if(neighborhood[inv_index] != land_grid_ptr[i*grid_ld[0]+j]) {
                                    updates.push_back({i, j, neighborhood[inv_index]});
                                }
                            }
                        }
                    }
                }
            }

            debug("After invasion %d", step);

            for(const auto& [i, j, val] : updates) {
                land_grid_ptr[i*grid_ld[0]+j] = val;
                bool unmask = true;
                for(int x = -1; x < 1; x++) {
                    for(int y = -1; y <= 1; y++) {
                        if((x != 0 || y != 0)) {
                            bool unmask = true;
                            for(int xx = -1; xx < 1; xx++) {
                                for(int yy = -1; yy <= 1; yy++) {
                                    if((xx != 0 || yy != 0)) {
                                        if(land_grid_ptr[(i+x+xx)*grid_ld[0]+j+y+yy] != val) {
                                            unmask = false;
                                            break;
                                        }
                                    }
                                }
                                if(!unmask) break;
                            }
                            land_mask_ptr[(i+x)*mask_ld[0]+j+y] = !unmask;
                        }
                    }
                }
            }

            debug("After updates %d", step);

            // renormalize every nstep steps
            if(step%endtime == 0) {
                cell_type land_grid_mean = 0;
                cell_type tmp_sum = 0;
                for(int i = lo[0]; i < hi[0]; i++) {
                    for(int j = lo[1]; j < hi[1]; j++) {
                        tmp_sum += land_grid_ptr[i*grid_ld[0]+j];
                    }
                    land_grid_mean += tmp_sum/(size*size);
                    tmp_sum = 0;
                }
                MPI_Allreduce(MPI_IN_PLACE, &land_grid_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for(int i = lo[0]; i < hi[0] && false; i++) {
                    for(int j = lo[1]; j < hi[1]; j++) {
                        land_grid_ptr[i*grid_ld[0]+j] /= land_grid_mean; // normalize grid
                    }
                }

                debug("After normalization %d", step);
            }

            debug("Completed step %d", step);

        }

        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);

        out_start_time = std::chrono::system_clock::now();
        FILE *fp;
        outfile = outfile + std::to_string(rep) + ".log";
        fp = fopen(outfile.c_str(), "w");
        GA_Print_file(fp, ga_land_grid);
        fclose(fp);
        out_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> out_time = out_end_time - out_start_time;
        println("Output run time = %fs", out_time.count());

    }

    end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;
    println("Total run time = %fs", total_time.count());
    println("Wrote results to file %s", outfile.c_str());
    fflush(stdout);

    GA_Destroy(ga_land_grid);
    GA_Destroy(ga_land_mask);

    GA_Terminate();
    MPI_Finalize();
}

