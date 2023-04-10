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
#include <iostream> // std::cout
#include <cxxopts.hpp> // to handle cmdline args
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

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

    cxxopts::Options options("sim", "Stocastic competition simulation");

    options.add_options()
        ("s,size", "grid size", cxxopts::value<int>()->default_value("500"))
        ("r,reps", "number of repetitions", cxxopts::value<int>()->default_value("1"))
        ("o,outfile", "output file name", cxxopts::value<std::string>()->default_value("out.csv"))
        ("h,help", "Print usage")
        ;

    auto result = options.parse(argc, argv);

    if(result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    outfile = result["outfile"].as<std::string>();
    size = result["size"].as<int>();
    nrep = result["reps"].as<int>();

    timescale = 100*(size*1.0/p);
    endtime = timescale/nsteps;

    start_time = std::chrono::system_clock::now();

    // grid for average across reps
    cell_type land_grid_data[size*size];
    cell_type* land_grid[size];
    bool land_mask_data[size*size];
    bool* land_mask[size];
    for(size_t i = 0; i < size; i++) {
        land_grid[i] = land_grid_data + i*size;
        land_mask[i] = land_mask_data + i*size;
    }

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

        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid[i][j] = 1;
                land_mask[i][j] = 0;
            }
        }

        // invasion rule variables
        cell_type neighborhood[8];
        cell_type inv[8];
        cell_type inv_sum = 0;
        int inv_index = 0;

        int i, j, row, col;
        cell_type min = 1;
        cell_type max = 1;

        for(int step = 0; step < timescale; step++) {

            // speciation rule
            int speciation_event_count = speciation_distribution(generator);
            std::set<int> spec_events;
            while(spec_events.size() < speciation_event_count) {
                int index = rand() % (size*size);
                if(!spec_events.count(index)) {
                    i = index % size;
                    j = index / size;
                    spec_events.insert(index);
                    int down = rand() % 2;
                    float ratio = (1+RanGen.Random()*mutsize);
                    if(down) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(RanGen.Random() <= probsuccess) {
                        land_grid[i][j] *= ratio;
                        if(land_grid[i][j] > max) {
                            max = land_grid[i][j];
                        }
                        if(land_grid[i][j] < min) {
                            min = land_grid[i][j];
                        }
                        for(int x = -1; x < 1; x++) {
                            for(int y = -1; y <= 1; y++) {
                                // wrap around
                                row = i+x;
                                col = j+y;
                                if(row >= size || row < 0)
                                    row = (row+size)%size;
                                if(col >= size || col < 0)
                                    col = (col+size)%size;
                                land_mask[row][col] = true;
                            }
                        }
                    }
                }
            }

            // invasion rule
            std::vector<cell_update> updates;
            for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++) {
                    if(land_mask[i][j]) {
                        double randval = RanGen.Random(); //random_float();
                        if(randval < invrate) {
                            inv_sum = 0;
                            inv_index = 0;
                            for(int x = -1; x < 1; x++) {
                                for(int y = -1; y <= 1; y++) {
                                    if((x != 0 || y != 0)) {
                                        row = i+x;
                                        col = j+y;
                                        if(row >= size || row < 0)
                                            row = (row+size)%size;
                                        if(col >= size || col < 0)
                                            col = (col+size)%size;
                                        neighborhood[inv_index] = land_grid[row][col];
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
                                    updates.push_back({i, j, neighborhood[inv_index]});
                                }
                            }
                        }
                    }
                }
            }

            for(const auto& [i, j, val] : updates) {
                land_grid[i][j] = val;
                bool unmask = true;
                for(int x = -1; x < 1; x++) {
                    for(int y = -1; y <= 1; y++) {
                        if((x != 0 || y != 0)) {
                            bool unmask = true;
                            for(int xx = -1; xx < 1; xx++) {
                                for(int yy = -1; yy <= 1; yy++) {
                                    row = i+x+xx;
                                    col = j+y+yy;
                                    if(row >= size || row < 0)
                                        row = (i+x+xx+size)%size;
                                    if(col >= size || col < 0)
                                        col = (j+y+yy+size)%size;
                                    if((xx != 0 || yy != 0)) {
                                        if(land_grid[row][col] != val) {
                                            unmask = false;
                                            break;
                                        }
                                    }
                                }
                                if(!unmask) break;
                            }
                            row = i+x;
                            col = j+y;
                            if(row >= size || row < 0)
                                row = (row+size)%size;
                            if(col >= size || col < 0)
                                col = (col+size)%size;
                            land_mask[row][col] = !unmask;
                        }
                    }
                }
            }

            // renormalize to prevent overflows
            if( std::numeric_limits<cell_type>::max() / (1+mutsize) < max ) {
                float land_grid_mean = 0;
                float tmp_sum = 0;
                for(int i = 0; i < size; i++) {
                    for(int j = 0; j < size; j++) {
                        tmp_sum += land_grid[i][j];
                    }
                    land_grid_mean += tmp_sum/(size*size);
                    tmp_sum = 0;
                }
                max = 1;
                min = 1;
                for(int i = 0; i < size; i++) {
                    for(int j = 0; j < size; j++) {
                        land_grid[i][j] /= land_grid_mean; // normalize grid
                        if(land_grid[i][j] > max) max = land_grid[i][j];
                        if(land_grid[i][j] < min) min = land_grid[i][j];
                    }
                }
                println("Normalized at step %d", step);
            }

        }

        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);

        out_start_time = std::chrono::system_clock::now();
        FILE *fp;
        outfile += std::to_string(rep) + ".log";
        fp = fopen(outfile.c_str(), "w");
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size-1; j++) {
                fprintf(fp, "%f", land_grid[i][j]);
                fprintf(fp, ", ");
            }
            fprintf(fp, "%f\n", land_grid[i][size-1]);
        }
        fclose(fp);
        out_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> out_time = out_end_time - out_start_time;
        println("Output run time = %fs", out_time.count());

    }

    end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;
    println("Total run time = %fs", total_time.count());
    println("Wrote reuslts to file %s", outfile.c_str());
    fflush(stdout);
}

