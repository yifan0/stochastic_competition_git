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
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;

#define println(...) { if(rank == 0) { printf(__VA_ARGS__); printf("\n"); } }
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
    MPI_Init(&argc, &argv);
    int nrep = 10;
    int size = 100;
    //int local_size, nblocks, size;
    double p = 0.1;
    double mutsize = 0.1;
    double specrate = 0.0001; // 1.0e-4
    double invrate = 0.2;
    int timescale = 100*size/p;
    int nsteps = 100;
    int endtime = timescale/nsteps;
    int rank, nprocs;
    std::string outfile = "output_";
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time, rep_start_time, rep_end_time, out_start_time, out_end_time;

    //MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
    outfile += std::to_string(rank) + ".log";
    timescale = 100*(size*1.0/p);
    endtime = timescale/nsteps;

    if(badargs) {
        println(">E Usage: %s\n", USAGE);
        MPI_Finalize();
        return 1;
    }

    int sqrt_nprocs = (int) std::sqrt(nprocs);
    if(sqrt_nprocs * sqrt_nprocs != nprocs) {
        if(rank == 0)
            fprintf(stderr, "Error: Number of MPI processes must be a perfect square.");
        MPI_Finalize();
        return 1;
    }

    println("Inputs:")
        println("\tMPI ranks = %d", nprocs)
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

    //print("Rank %d\n", rank);


    // Calculate the size of local submatrix
    int sub_size = size / sqrt_nprocs;
    int sub_size_remainder = size % sqrt_nprocs;
    int sub_size_start = (rank % sqrt_nprocs) * sub_size;
    if(rank % sqrt_nprocs > sub_size_remainder) {
        sub_size_start += sub_size_remainder;
        sub_size -= 1;
    }
    int sub_size_end = sub_size_start + sub_size - 1;
    if(rank % sqrt_nprocs == sub_size_remainder - 1) {
        sub_size += 1;
    }

    println("Size for each process = %d", sub_size);

    start_time = std::chrono::system_clock::now();

    // Random number generation
    CRandomSFMT1 RanGen((int)time(NULL)); // Agner Combined generator
    // distribution for speciation events
    std::default_random_engine generator;
    std::binomial_distribution<int> speciation_distribution(sub_size*sub_size, specrate);

    // grid for land data
    boost::numeric::ublas::matrix<cell_type> land_grid(sub_size+2, sub_size+2);

    MPI_Win win;
    MPI_Win_create(&land_grid(0,0), (sub_size+2) * (sub_size+2) * sizeof(cell_type), sizeof(cell_type), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    MPI_Win_fence(0, win);

    // Exchange ghost cells with neighboring processes using one-sided communication
    int top_proc = rank - sqrt_nprocs;
    int bottom_proc = rank + sqrt_nprocs;
    int left_proc = rank - 1;
    int right_proc = rank + 1;
    if (rank % sqrt_nprocs == 0) {
        left_proc = MPI_PROC_NULL;
    }
    if (rank % sqrt_nprocs == sqrt_nprocs - 1) {
        right_proc = MPI_PROC_NULL;
    }
    if (rank < sqrt_nprocs || nprocs == 1) {
        top_proc = MPI_PROC_NULL;
    }
    if (rank >= nprocs - sqrt_nprocs) {
        bottom_proc = MPI_PROC_NULL;
    }

    for(int rep = 0; rep < nrep; rep++) {
        rep_start_time = std::chrono::system_clock::now();

        // set all of land_grid_data to 1
        for(size_t i = 0; i < sub_size+2; i++) {
            for(size_t j = 0; j < sub_size+2; j++) {
                land_grid(i,j) = 1;
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
                int index = rand() % (sub_size*sub_size);
                if(!spec_events.count(index)) {
                    int i = index % sub_size;
                    int j = index / sub_size;
                    spec_events.insert(index);
                    int down = rand() % 2;
                    float ratio = (1+RanGen.Random()*mutsize);
                    if(down) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(RanGen.Random() <= probsuccess) {
                        i++; j++; // to account for ghost cells
                        land_grid(i,j) *= ratio;
                        if(i == 0 && top_proc != MPI_PROC_NULL) {
                            MPI_Put(&land_grid(i,j), 1, MPI_DOUBLE, top_proc, (sub_size+2)*(sub_size+1)+j, 1, MPI_DOUBLE, win);
                        }
                        if(i == sub_size-1 && bottom_proc != MPI_PROC_NULL) {
                            MPI_Put(&land_grid(i,j), 1, MPI_DOUBLE, bottom_proc, j, 1, MPI_DOUBLE, win);
                        }
                        if(j == 0 && left_proc != MPI_PROC_NULL) {
                            MPI_Put(&land_grid(i,j), 1, MPI_DOUBLE, left_proc, (sub_size+2)*(i+1)-1, 1, MPI_DOUBLE, win);
                        }
                        if(j == sub_size-1 && right_proc != MPI_PROC_NULL) {
                            MPI_Put(&land_grid(i,j), 1, MPI_DOUBLE, right_proc, (sub_size+2)*(i), 1, MPI_DOUBLE, win);
                        }
                    }
                }
            }

            MPI_Win_fence(0, win);

            // invasion rule
            std::vector<cell_update> updates;
            for(int i = 1; i < sub_size+1; i++) {
                for(int j = 1; j < sub_size+1; j++) {
                    double randval = RanGen.Random(); //random_float();
                    if(randval < invrate) {
                        inv_sum = 0;
                        inv_index = 0;
                        for(int x = -1; x <= 1; x++) {
                            for(int y = -1; y <= 1; y++) {
                                if((x != 0 || y != 0) && (i+x >= 1 || left_proc != MPI_PROC_NULL) && (i+x < sub_size+1 || right_proc != MPI_PROC_NULL) && (j+y >= 1 || top_proc != MPI_PROC_NULL) && (j+y < sub_size+1 || bottom_proc != MPI_PROC_NULL)) {
                                    neighborhood[inv_index] = land_grid(i+x,j+y);
                                    inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+land_grid(i,j)*(1-p));
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
                            if(neighborhood[inv_index] != land_grid(i,j)) {
                                //land_grid[i][j] = neighborhood[inv_index];
                                updates.push_back({i, j, neighborhood[inv_index]});
                            }
                        }
                    }
                }
            }
            for(const auto& [i, j, val] : updates) {
                land_grid(i,j) = val;
                if(i == 0 && top_proc != MPI_PROC_NULL) {
                    MPI_Put(&land_grid(i,j), 1, MPI_DOUBLE, top_proc, (sub_size+2)*(sub_size+1)+j, 1, MPI_DOUBLE, win);
                }
                if(i == sub_size-1 && bottom_proc != MPI_PROC_NULL) {
                    MPI_Put(&land_grid(i,j), 1, MPI_DOUBLE, bottom_proc, j, 1, MPI_DOUBLE, win);
                }
                if(j == 0 && left_proc != MPI_PROC_NULL) {
                    MPI_Put(&land_grid(i,j), 1, MPI_DOUBLE, left_proc, (sub_size+2)*(i+1)-1, 1, MPI_DOUBLE, win);
                }
                if(j == sub_size-1 && right_proc != MPI_PROC_NULL) {
                    MPI_Put(&land_grid(i,j), 1, MPI_DOUBLE, right_proc, (sub_size+2)*(i), 1, MPI_DOUBLE, win);
                }
            }
            
            //MPI_Win_fence(0, win);

            // renormalize every nstep steps
            // normalize the simulated result so that the fittest species will not have a very high effective population so that it always outcompetes other species
            if(step%nsteps == 0) {
                // Calculate a global average. This is not precise because it double-counts the boundaries, but it's a rough estimate for normalization
                cell_type local_average = 0;
                for(size_t i = 1; i < sub_size+1; i++) {
                    cell_type row_sum = 0;
                    for(size_t j = 1; j < sub_size+1; j++) {
                        row_sum += land_grid(i,j);
                    }
                    local_average += row_sum/(sub_size*sub_size);
                }
                //cell_type sum = std::accumulate(land_grid.begin1(), land_grid.end1(), 0.0);
                //cell_type local_average = sum / (land_grid.size1()*land_grid.size2());
                cell_type global_average = 0;
                MPI_Allreduce(&local_average, &global_average, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                global_average /= nprocs;

                // Normalize based on the global average
                for(size_t i = 0; i < sub_size+2; i++) {
                    for(size_t j = 0; j < sub_size+2; j++) {
                        land_grid(i,j) /= global_average;
                    }
                }
                //std::transform(land_grid.begin1(), land_grid.end1(), land_grid.begin1(), [global_average](double v) { return v / global_average; });
            }

        }

        // Normalize
        cell_type local_average = 0;
        for(size_t i = 1; i < sub_size+1; i++) {
            cell_type row_sum = 0;
            for(size_t j = 1; j < sub_size+1; j++) {
                row_sum += land_grid(i,j);
            }
            local_average += row_sum/(sub_size*sub_size);
        }
        cell_type global_average = 0;
        MPI_Allreduce(&local_average, &global_average, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_average /= nprocs;
        for(size_t i = 0; i < sub_size+2; i++) {
            for(size_t j = 0; j < sub_size+2; j++) {
                land_grid(i,j) /= global_average;
            }
        }

        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);

        // Output for this rep
        out_start_time = std::chrono::system_clock::now();
        if(argc > 3)
            outfile = argv[3];
        else
            outfile = "output_";
        outfile += "rep" + std::to_string(rep) + "_rank" + std::to_string(rank) + ".log";
        FILE *fp;
        fp = fopen(outfile.c_str(), "w");
        print("Writing to file %s\n", outfile.c_str());
        for(int i = 1; i < sub_size+1; i++) {
            for(int j = 1; j < sub_size; j++) {
                fprintf(fp, "%f", land_grid(i,j));
                fprintf(fp, ", ");
            }
            fprintf(fp, "%f\n", land_grid(i,sub_size));
        }
        fflush(fp);
        fclose(fp);
        out_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> out_time = out_end_time - out_start_time;
        println("Output run time = %fs", out_time.count());
        fflush(stdout);


    }

    end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;
    println("Total run time = %fs", total_time.count());
    fflush(stdout);

    MPI_Finalize();
    return 0;
}

