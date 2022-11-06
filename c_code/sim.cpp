#include <cstdio>
#include <vector>
#include <random>
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

enum DIRECTION { NW=0, N=1, NE=2, W=10, HERE=11, E=12, SW=20, S=21, SE=22 };

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

int main(int argc, char* argv[]){
    int nrep = 1;
    int size = 128;
    double p = 0.1;
    double mutsize = 0.1;
    double specrate = 0.0001; // 1.0e-4
    int timescale = 100/p*size;
    int nsteps = 100;
    int endtime = timescale/nsteps;
    double land_grid[size][size];

    println("Inputs:")
    println("\trepetitions = %d", nrep)
    println("\tsize = %d%s%d", size, "x", size)
    println("\tindividuals per patch = %f", 1/p)
    println("\tmutation size = %f", mutsize)
    println("\tspeciation rate = %f", specrate)
    println("\ttimescale = %d", timescale)
    println("\tnsteps = %d", nsteps)
    println("\tend time = %d", endtime)

    for(int rep = 0; rep < nrep; rep++) {
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
                                if(x != 0 || y != 0) {
                                    neighborhood[inv_index] = land_grid[i+x][j+y];
                                    inv[inv_index] = p*neighborhood[inv_index]/(p*neighborhood[inv_index]+land_grid[i][j]*(1-p));
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
                            // TODO: this should be a neighbor value
                            land_grid[i][j] = neighborhood[inv_index];
                        }
                    }
                }
            }
        }
        // normalize the simulated result so that the fittest species will not have a very high effective population so that it always outcompetes other species
        float land_grid_mean = 0;
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid_mean += land_grid[i][j]/(size*size);
            }
        }
        println("Grid mean = %f", land_grid_mean);
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid[i][j] /= land_grid_mean;
            }
        }
        // TODO: write to file
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                print("%f", land_grid[i][j]);
                print(" ");
            }
            println("");
        }
    }
    // TODO: average all of the grids

}

