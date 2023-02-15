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
#include <bcl/bcl.hpp>
#include <bcl/containers/DMatrix.hpp>
using namespace std;

#define USAGE \
"bclsim gridsize nreps outfile"

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

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;
float random_float(float max=1.0) {
    // returns a random value between 0 and 1
    return (FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN))))*max;
}



int main(int argc, char* argv[]){
    BCL::init();
    size_t nrep = 10;
    size_t size = 100;
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
    if(argc > 4) badargs = true;
    if(argc > 1 && sscanf(argv[1],"%ld",&size) != 1)
        badargs = true;
    if(size < 1) badargs = true;
    if (argc > 2 && sscanf(argv[2],"%ld",&nrep) != 1)
        badargs = true;
    if(nrep < 1) badargs = true;
    if(argc > 3)
       outfile = argv[3];
    timescale = 100*(size*1.0/p);
    endtime = timescale/nsteps;

    if(badargs) {
        BCL::print(">E Usage: %s\n", USAGE);
        BCL::finalize();
        return 1;
    }
    
    start_time = std::chrono::system_clock::now();

    // grid for average across reps
    BCL::DMatrix<float> mean_grid({size, size});
    BCL::DMatrix<float> rep_grid({size, size});

    BCL::print("Inputs:\n");
    BCL::print("\trepetitions = %d\n", nrep);
    BCL::print("\tsize = %d%s%d\n", size, "x", size);
    BCL::print("\tindividuals per patch = %f\n", 1/p);
    BCL::print("\tmutation size = %f\n", mutsize);
    BCL::print("\tspeciation rate = %f\n", specrate);
    BCL::print("\tinvasion rate = %f\n", invrate);
    BCL::print("\ttimescale = %d\n", timescale);
    BCL::print("\tnsteps = %d\n", nsteps);
    BCL::print("\tend time = %d\n", endtime);
    BCL::print("\n");

    // Random number generation
    CRandomSFMT1 RanGen((int)time(NULL)); // Agner Combined generator
    // distribution for speciation events
    std::default_random_engine generator;
    std::binomial_distribution<int> speciation_distribution(size*size, specrate);

    // shape() gives the total size of the matrix
    // tile_shape() gives the dimensions of the block assigned to each process
    // grid_shape() gives the number of processes by the number of processes
    // For the sizes of the edge blocks, need to get info from tile_rank() and then tile_shape({i,j})
    BCL::print("pm() = %d, pn() = %d\n", rep_grid.pm(), rep_grid.pn());

    size_t local_rows = rep_grid.tile_shape()[0];
    size_t local_cols = rep_grid.tile_shape()[1];
    printf("Rank %d: rows = %d, cols = %d. Size = %d\n", BCL::rank(), local_rows, local_cols, rep_grid.tile_size());

    if(BCL::rank() == 0)
    rep_grid.print_info();

    // invasion rule variables
    cell_type neighborhood[8];
    cell_type inv[8];
    cell_type inv_sum = 0;
    int inv_index = 0;

    for(int rep = 0; rep < nrep; rep++) {
        rep_start_time = std::chrono::system_clock::now();

        rep_grid = 1;

        for(int step = 0; step < timescale; step++) {
            int speciation_event_count = speciation_distribution(generator);
        }
        /*
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
        std::map<int, cell_type> invasion_event_archive;
        for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
        double randval = RanGen.Random(); //random_float();
        if(randval < invrate) {
        inv_sum = 0;
        inv_index = 0;
        for(int x = -1; x < 1; x++) {
        for(int y = -1; y <= 1; y++) {
        if((x != 0 || y != 0) && i+x >= 0 && i+x < size && j+y >= 0 && j+y < size) {
        if(invasion_event_archive.find( (i+x)*size+j+y ) != invasion_event_archive.end() )
        neighborhood[inv_index] = invasion_event_archive[(i+x)*size+j+y];
        else
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
        land_grid[i][j] = neighborhood[inv_index];
        invasion_event_archive[i*size + j] = neighborhood[inv_index];
        }
    }
    }
    }
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
    */
        rep_end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
    BCL::print("Rep %d run time = %fs\n", rep, rep_time.count());
    }

    out_start_time = std::chrono::system_clock::now();
    /*
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
       */
    out_end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> out_time = out_end_time - out_start_time;
    BCL::print("Output run time = %fs\n", out_time.count());
    end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;
    BCL::print("Total run time = %fs\n", total_time.count());
    BCL::print("Wrote results to file %s\n", outfile.c_str());

    BCL::finalize();
}

