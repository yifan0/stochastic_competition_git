#include <cstdio>
#include <random> // rand()
#include <string.h> // memcpy()
#include <queue> // queue
#include <omp.h>
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

enum DIRECTION { NW=0, N=1, NE=2, W=10, HERE=11, E=12, SW=20, S=21, SE=22 };

typedef double cell_type;

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;
float random_float(float max=1.0) {
    // returns a random value between 0 and 1
    return (FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN))))*max;
}

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
    int size = 512;
    double p = 0.1;
    double mutsize = 0.1;
    double specrate = 0.0001; // 1.0e-4
    int timescale = 100*size/p;
    int nsteps = 100;
    int endtime = timescale/nsteps;
    std::string outfile = "cpptest128grid_result.txt";

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
    println("\ttimescale = %d", timescale)
    println("\tnsteps = %d", nsteps)
    println("\tend time = %d", endtime);
    println("");

    for(int rep = 0; rep < nrep; rep++) {
        //println("Number of threads = %d", omp_get_num_threads());
        // grid for land data
        // strategy from: https://stackoverflow.com/questions/46354262/how-to-dynamically-allocate-a-contiguous-2d-array-in-c
        std::vector<cell_type> land_grid_data(size*size);
        std::vector<cell_type*> land_grid_arrays;
        for(int i = 0; i != size*size; i += size) {
            land_grid_arrays.push_back(land_grid_data.data() + i);
        }
        cell_type** land_grid = land_grid_arrays.data();

        // grid for old land data to save from previous iteration
        std::vector<cell_type> old_land_grid_data(size*size);
        std::vector<cell_type*> old_land_grid_arrays;
        for(int i = 0; i != size*size; i += size) {
            old_land_grid_arrays.push_back(old_land_grid_data.data() + i);
        }
        cell_type** old_land_grid = old_land_grid_arrays.data();
        
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid[i][j] = 1;
            }
        }
        
        for(int step = 0; step < endtime; step++) {
            // do the speciation rule
            for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++) {
                    double randval = random_float();
                    if(randval < specrate) {
                        int down = rand() % 2; // if false, increase effective population, otherwise decrease
                        float ratio = (1+random_float()*mutsize);
                        if(down) {
                            ratio = 1/ratio;
                        }
                        float probsuccess = p*ratio/(p*(ratio-1)+1);
                        if(randval <= specrate*probsuccess) {
                            land_grid[i][j] *= ratio;
                        }
                    }
                }
            }
            // do the neighbors rule
           
            // save updated local values as old_land_grid and fill land_grid with new vals
            swap(land_grid, old_land_grid);
            
            for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++) {
                    double randval = random_float();
                    if (randval < 0.2) {
                        // find neighbor values
                        float neighborhood[8];
                        float inv[8];
                        float inv_sum = 0;
                        int inv_index = 0;
                        for(int x = -1; x <= 1; x++) {
                            for(int y = -1; y <= 1; y++) {
                                if((x != 0 || y != 0) && i+x >= 0 && i+x < size && j+y >= 0 && j+y < size) {
                                    //neighborhood[inv_index] = land_grid[i+x][j+y];
                                    neighborhood[inv_index] = old_land_grid[i+x][j+y];
                                    //inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+land_grid[i][j]*(1-p));
                                    inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+old_land_grid[i][j]*(1-p));
                                    inv_sum += inv[inv_index];
                                    inv_index++;
                                }
                            }
                        }
                        if(randval <= inv_sum/8) {
                            // weight the random number
                            float ran = random_float(inv_sum);
                            inv_index = 0;
                            while(ran > inv[inv_index]) {
                                ran -= inv[inv_index];
                                inv_index++;
                            }
                            //cellUpdates.push(newCellUpdateRecord(i, j, neighborhood[inv_index]));
                            land_grid[i][j] = neighborhood[inv_index];
                        } else {
                            land_grid[i][j] = old_land_grid[i][j];
                        }
                    } else {
                        land_grid[i][j] = old_land_grid[i][j];
                    }
                }
            }

            /*
            cellUpdateRecord rec;
            while(!cellUpdates.empty()) {
                rec = cellUpdates.front();
                land_grid[rec.y][rec.x] = rec.val;            
                cellUpdates.pop();
            }*/
        }

        // normalize the simulated result so that the fittest species will not have a very high effective population so that it always outcompetes other species
        float land_grid_mean = 0;
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid_mean += land_grid[i][j]/(size*size);
            }
        }
        println("Grid mean = %f for rep %d", land_grid_mean, rep);
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid[i][j] /= land_grid_mean; // normalize grid
                mean_grid[i][j] += land_grid[i][j]/nrep; // add contribution to average across reps
            }
        }
    }

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
    println("Wrote reuslts to file %s", outfile.c_str());
}

