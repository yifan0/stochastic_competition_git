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
#include <stdio.h> //remove file
#include "width1D.h"
using namespace std;

#define println(...) { if(me == 0) { printf(__VA_ARGS__); printf("\n"); } }
#define debug(...) { if(true) { printf("Rank %d: ", me); printf(__VA_ARGS__); printf("\n"); GA_Sync(); } }
#define print(...) { printf(__VA_ARGS__); }

#define   NDIM         2
#define   GHOSTS       2

typedef double cell_type;
typedef std::tuple<int, int, cell_type> cell_update;

int main(int argc, char* argv[]){
    ofstream myfile;
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
    string width_file;
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time, rep_start_time, rep_end_time, out_start_time, out_end_time, step_start_time, step_end_time, save_start_time, save_end_time, width_start_time, width_end_time;
    int me, nprocs;
    int dims[NDIM];
    int grid_ld[NDIM];
    int lo[NDIM], hi[NDIM];
    int ghost_grid_ld[NDIM-1];
    int ghost_dims[NDIM];
    int ghost_width[NDIM];
    double* land_grid_ptr, *ghost_grid_ptr;
    int heap=3000000, stack=3000000; // TODO: find out if these are good values
    FILE *fp;

    MPI_Init(&argc, &argv);
    GA_Initialize();
    char init_err_msg[] = "MA_init failed";
    if(! MA_init(C_DBL, stack, heap) ) GA_Error(init_err_msg,stack+heap);

    me = GA_Nodeid();
    nprocs = GA_Nnodes();

    cxxopts::Options options("mpirun -n N ./ga_sim", "Stocastic competition simulation with GlobalArray parallelism");

    options.add_options()
        ("s,size", "grid size", cxxopts::value<int>()->default_value("500"))
        ("c,specrate", "speciation rate", cxxopts::value<double>()->default_value("0.0001"))
        ("r,reps", "number of repetitions", cxxopts::value<int>()->default_value("1"))
        ("m,mutsize", "maximum change in mutation event", cxxopts::value<double>()->default_value("0.1"))
        ("o,outfile", "output file prefix", cxxopts::value<std::string>()->default_value("out"))
        ("a,width_file", "save file prefix", cxxopts::value<std::string>()->default_value("out"))
        ("h,help", "Print usage")
        ;

    auto result = options.parse(argc, argv);

    if(result.count("help")) {
        if(me == 0)
            std::cout << options.help() << std::endl;
        MPI_Finalize();
        exit(0);
    }
    size = result["size"].as<int>();
    nrep = result["reps"].as<int>();
    specrate = result["specrate"].as<double>();
    mutsize = result["mutsize"].as<double>();
    width_file = result["width_file"].as<string>();
    outfile = result["outfile"].as<string>();
    nsteps = 10;
    timescale = 10*(((int)sqrt(size))*1.0/p);
    endtime = timescale/nsteps;
    double width_step[nsteps][23];

    // check that nprocs is a square to prevent errors
    int sqrt_nprocs = sqrt(nprocs);
    if(sqrt_nprocs*sqrt_nprocs != nprocs) {
        println("Number of processes must form a square");
        MPI_Finalize();
        exit(0);
    }

    start_time = std::chrono::system_clock::now();

    // grid for average across reps
    GA_Mask_sync(0, 0); // turns off sync when updating ghosts
    dims[0] = 1;
    dims[1] = size;
    ghost_width[0] = 0;
    ghost_width[1] = 2;
    char land_grid_name[] = "land grid";
    int ga_land_grid = NGA_Create_ghosts(C_DBL, NDIM, dims, ghost_width, land_grid_name, NULL);

    NGA_Distribution(ga_land_grid, me, lo, hi);
    NGA_Access(ga_land_grid, lo, hi, &land_grid_ptr, grid_ld);

    println("Inputs:")
        println("\trepetitions = %d", nrep)
        println("\tsize = %d%s%d", dims[0], "x", dims[1])
        println("\tindividuals per patch = %f", 1/p)
        println("\tmutation size = %f", mutsize)
        println("\tspeciation rate = %.2e", specrate)
        println("\tinvasion rate = %f", invrate)
        println("\ttimescale = %d", timescale)
        println("\tnsteps = %d", nsteps)
        println("\tend time = %d", endtime);
        println("\toutfile = %s", outfile.c_str());
        println("\twidth_file = %s", width_file.c_str());
    println("");
    fflush(stdout);

    // Random number generation
    CRandomSFMT1 RanGen(10); // Agner Combined generator
    int local_cols = (hi[0]-lo[0]+1);
    int local_rows = (hi[1]-lo[1]+1);
    println("local_cols = %d", local_cols);
    println("local_rows = %d", local_rows);
    int local_area = local_cols*local_rows;
    // double land_mask_data[local_area];
    // double* land_mask[local_rows];
    // for(size_t i = 0; i < local_rows; i++) {
    //     land_mask[i] = land_mask_data + i*local_cols;
    // }
    // distribution for speciation events
    std::default_random_engine generator;
    std::binomial_distribution<int> speciation_distribution(local_area, specrate);

    for(int rep = 0; rep < nrep; rep++) {
        rep_start_time = std::chrono::system_clock::now();
        step_start_time = std::chrono::system_clock::now();
        double save_total = 0;
        double step_total = 0;
        double width_total = 0;
        int zero = 0;
        double one = 1;
        GA_Fill(ga_land_grid, &one);
        GA_Update_ghosts(ga_land_grid);
        NGA_Access_ghosts(ga_land_grid, ghost_dims, &ghost_grid_ptr, ghost_grid_ld);
        println("ghost_grid_ld = %d", ghost_grid_ld[0]);
        // zero out mask except local edges
        // memset(land_mask_data, 0, sizeof(*land_mask_data));
        // for(size_t i = 0; i < local_rows; i++) {
        //     land_mask[i][0] = 1;
        //     land_mask[i][local_cols-1] = 1;
        // }
        // for(size_t i = 0; i < local_cols; i++) {
        //     land_mask[0][i] = 1;
        //     land_mask[local_rows-1][i] = 1;
        // }
        double land_mask[local_rows];
        for (int i=0; i<local_rows; i++){
            land_mask[i] = 0;
        }
        land_mask[0] = 1;
        land_mask[local_rows-1] = 1;
        int boundary_count = 0;
        for (int i=0; i<local_rows; i++){
            if (land_mask[i] != 0){
                boundary_count += 1;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &boundary_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        println("number of cell on boundary in the beginning = %d", boundary_count);
        // invasion rule variables
        cell_type neighborhood[2];
        cell_type inv[2];
        cell_type inv_sum = 0;
        int inv_index = 0;
        int row, col;

        for(int step = 0; step < timescale; step++) {

            // speciation rule
            int speciation_event_count = speciation_distribution(generator);
            std::set<int> spec_events;
            while(spec_events.size() < speciation_event_count) {
                int index = rand() % ((hi[0]-lo[0]+1)*(hi[1]-lo[1]+1));
                if(!spec_events.count(index)) {
                    int i = index / local_rows;
                    int j = index % local_rows;
                    spec_events.insert(index);
                    int down = rand() % 2;
                    float ratio = (1+RanGen.Random()*mutsize);
                    if(down) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(RanGen.Random() <= probsuccess) {
                        ghost_grid_ptr[(i+ghost_width[0])*grid_ld[0]+j+ghost_width[1]] *= ratio;
                    } else {
                        continue;
                    }

                    // set mask for [i] and surrounding cells [i-1], [i+1] unless on edge (but only if the speciation is fixed)
                    for(int y = -1; y <= 1; y++) {
                        row = i;
                        col = j+y;
                        // if(row > local_cols-1 || row < 0)
                            // continue;
                        if(col >= local_rows-1 || col <= 0)
                            continue;
                        cell_type local_max = ghost_grid_ptr[(row+ghost_width[0])*grid_ld[0]+col+ghost_width[1]];
                        for(int yy = -1; yy <= 1; yy++) {
                            int neighbor_row = row;
                            int neighbor_col = col+yy;
                            local_max = std::max(local_max, ghost_grid_ptr[(neighbor_row+ghost_width[0])*grid_ld[0]+neighbor_col+ghost_width[1]]);
                        }
                        land_mask[col] = local_max*p/(local_max*p+ghost_grid_ptr[(row+ghost_width[0])*grid_ld[0]+col+ghost_width[1]]*(1-p));
                    }
                }
            }

            if(step%10 == 0)
                GA_Update_ghosts(ga_land_grid);

            // invasion rule
            std::vector<cell_update> updates;
            for(int i = 0; i < local_cols; i++) {
                for(int j = 0; j < local_rows; j++) {
                    if(land_mask[j] != 0) {
                        double randval = RanGen.Random(); //random_float();
                        if(randval < land_mask[j]) {
                            inv_sum = 0;
                            inv_index = 0;
                            for(int y = -1; y <= 1; y++) {
                                if(y != 0) {
                                    neighborhood[inv_index] = ghost_grid_ptr[(i+ghost_width[0])*grid_ld[0]+j+ghost_width[1]+y];
                                    inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+ghost_grid_ptr[(i+ghost_width[0])*grid_ld[0]+j+ghost_width[1]]*(1-p));
                                    inv_sum += inv[inv_index];
                                    inv_index++;
                                }
                            }

                            if(randval <= inv_sum/2) {
                                // Get random element with weighted probabilities
                                double weighted_rand = RanGen.Random()*inv_sum;
                                inv_index = 0;
                                while(weighted_rand > inv[inv_index]) {
                                    weighted_rand -= inv[inv_index];
                                    inv_index++;
                                }
                                if(neighborhood[inv_index] != ghost_grid_ptr[(i+ghost_width[0])*grid_ld[0]+j+ghost_width[1]]) {
                                    updates.push_back({i, j, neighborhood[inv_index]});
                                }
                            }
                        }
                    }
                }
            }

            for(const auto& [i, j, val] : updates) {
                ghost_grid_ptr[(i+ghost_width[0])*grid_ld[0]+j+ghost_width[1]] = val;
                for(int y = -1; y <= 1; y++) {
                    row = i;
                    col = j+y;
                    // if(row >= local_cols-1 || row <= 0)
                        // continue;
                    if(col >= local_rows-1 || col <= 0)
                        continue;
                    cell_type local_max = 0;
                    for(int yy = -1; yy <= 1; yy++) {
                        row = i;
                        col = j+y+yy;
                        // if(row >= local_cols-1 || row <= 0)
                            // continue;
                        if(col >= local_rows-1 || col <= 0)
                            continue;
                        if(ghost_grid_ptr[(row+ghost_width[0])*grid_ld[0]+col+ghost_width[1]] != val && ghost_grid_ptr[(row+ghost_width[0])*grid_ld[0]+col+ghost_width[1]] > local_max) {
                            local_max = ghost_grid_ptr[(row+ghost_width[0])*grid_ld[0]+col+ghost_width[1]];
                        }
                    }
                    row = i;
                    col = j+y;
                    if(local_max == 0) {
                        land_mask[col] = 0;
                    }
                    else {
                        land_mask[col] = local_max*p/(local_max*p+ghost_grid_ptr[(row+ghost_width[0])*grid_ld[0]+col+ghost_width[1]]*(1-p));
                    }
                }
            }

            // renormalize every nstep steps
            if(step%endtime == endtime-1) {
                println("step = %d", step);
                cell_type land_grid_mean = 0;
                for(int i = 0; i <= hi[0]-lo[0]; i++) {
                    for(int j = 0; j <= hi[1]-lo[1]; j++) {
                        land_grid_mean += ghost_grid_ptr[(i+ghost_width[0])*grid_ld[0]+j+ghost_width[1]];
                    }
                }
                land_grid_mean /= (dims[0]*dims[1]);
                MPI_Allreduce(MPI_IN_PLACE, &land_grid_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for(int i = 0; i <= hi[0]-lo[0]; i++) {
                    for(int j = 0; j <= hi[1]-lo[1]; j++) {
                        ghost_grid_ptr[(i+ghost_width[0])*grid_ld[0]+j+ghost_width[1]] /= land_grid_mean; // normalize grid
                    }
                }
                GA_Update_ghosts(ga_land_grid);
                outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + "_checkpoint" + std::to_string((int)floor(step/endtime)) + ".csv";
                save_start_time = std::chrono::system_clock::now();
                fp = fopen(outfile.c_str(), "w");
                GA_Print_csv_file(fp, ga_land_grid);
                fclose(fp);
                save_end_time = std::chrono::system_clock::now();
                std::chrono::duration<double> save_time = save_end_time - save_start_time;
                save_total += save_time.count();
                println("save time = %fs", save_time.count());
                println("Output with progress %d% (step = %d) to file %s", step/endtime, step, outfile.c_str());
                width_start_time = std::chrono::system_clock::now();
                // if (me == 0){
                    // double *pwidth = width1D_simple(outfile.c_str());
                    double *pwidth = width1D_simple_parallel(outfile.c_str());
                    for (int i=0; i<23; i++){
                        width_step[(int) step/endtime][i] = *(pwidth+i);
                    }
                // }
                MPI_Barrier(MPI_COMM_WORLD);
                width_end_time = std::chrono::system_clock::now();
                std::chrono::duration<double> width_time = width_end_time - width_start_time;
                width_total += width_time.count();
                println("width time = %fs", width_time.count());
                step_end_time = std::chrono::system_clock::now();
                std::chrono::duration<double> step_time = step_end_time - step_start_time;
                step_total += step_time.count();
                println("step run time = %fs", step_time.count());
                step_start_time = std::chrono::system_clock::now();
                MPI_Barrier(MPI_COMM_WORLD);
                remove(outfile.c_str());
            }
        }
        println("avg step run time = %fs", step_total/nsteps);
        println("avg save run time = %fs", save_total/nsteps);
        println("avg width run time = %fs", width_total/nsteps);
        println("avg step run time pure = %fs", (step_total-save_total-width_total)/nsteps);
        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);

        out_start_time = std::chrono::system_clock::now();
        outfile = result["outfile"].as<std::string>() + "_rep" + std::to_string(rep) + ".csv";
        fp = fopen(outfile.c_str(), "w");
        GA_Print_csv_file(fp, ga_land_grid);
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

    GA_Terminate();
    MPI_Finalize();
    myfile.open (width_file.c_str());
    for (int i=0; i<nsteps; i++){
        for (int j=0; j<20; j++){
            myfile << width_step[i][j] << ",";
        }
        myfile << endl;
    }
    myfile.close();
    println("Saved to %s", width_file.c_str());
}

