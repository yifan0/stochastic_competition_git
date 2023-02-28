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

struct speciation_tree_node {
    cell_type val;
    size_t i;
    size_t j;
    size_t time;
    vector<speciation_tree_node*> children;
    speciation_tree_node(cell_type v, size_t t, size_t x, size_t y) : val(v), time(t), i(x), j(y) {
    }
};

typedef std::tuple<int, int, speciation_tree_node*> cell_update;

speciation_tree_node* find_node(speciation_tree_node* root, cell_type target) {
    if(root->val == target) {
        return root;
    }

    for(auto child : root->children) {
        speciation_tree_node* result = find_node(child, target);
        if(result != nullptr) return result;
    }

    return nullptr;
}

// caller is responsible for adding a semicolon at the end
string toString(speciation_tree_node* root) {
    // for a leaf, just return the label
    if(root->children.size() == 0) {
        return std::to_string(root->val);
    }

    // for a parent, group the children in () and add label
    string output = "(";
    for(int i = 0; i < root->children.size(); i++) {
        output += toString(root->children[i]);
        if(i < root->children.size()-1) output += ",";
    }
    output += ")";
    output += std::to_string(root->val);
    return output;
}

void delete_tree(speciation_tree_node* root) {
    if(root == nullptr) return;
    for(auto child : root->children) {
        delete_tree(child);
    }
    delete root;
}

size_t get_depth(speciation_tree_node* root) {
    if(root == nullptr) return 0;
    size_t depth = 0;
    for(auto child : root->children) {
        depth = max(depth, get_depth(child));
    }
    return depth+1;
}

int main(int argc, char* argv[]){
    int nrep = 10;
    int size = 100;
    double p = 0.1;
    double mutsize = 0.1;
    double specrate = 0.0001; // 1.0e-4
    double invrate = 0.2;
    invrate = 0;
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
    vector<speciation_tree_node*> land_grid_data(size*size);
    vector<speciation_tree_node**> land_grid(size);
    //cell_type land_grid_data[size*size];
    //cell_type* land_grid[size];
    bool land_mask_data[size*size];
    bool* land_mask[size];
    for(size_t i = 0; i < size; i++) {
        land_grid[i] = &land_grid_data[i*size];
        land_mask[i] = &land_mask_data[i*size];
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

        // speciation rule variables
        speciation_tree_node* speciation_root = new speciation_tree_node(1, 0, 0, 0); // root has value 1, the starting value

        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                land_grid[i][j] = speciation_root;
                land_mask[i][j] = 0;
            }
        }

        // invasion rule variables
        speciation_tree_node* neighborhood[8];
        cell_type inv[8];
        cell_type inv_sum = 0;
        int inv_index = 0;

        for(int step = 0; step < timescale; step++) {

            // speciation rule
            int speciation_event_count = speciation_distribution(generator);
            std::set<int> spec_events;
            while(spec_events.size() < speciation_event_count) {
                int index = RanGen.IRandom(0,size*size-1);
                if(!spec_events.count(index)) {
                    int i = index % size;
                    int j = index / size;
                    spec_events.insert(index);
                    double ratio = 1+RanGen.Random()*mutsize;
                    if(RanGen.IRandom(0,1)) {
                        ratio = 1/ratio;
                    }
                    float probsuccess = p*ratio/(p*(ratio-1)+1);
                    if(RanGen.Random() <= probsuccess) {
                        speciation_tree_node* parent = land_grid[i][j];

                        land_grid[i][j] = new speciation_tree_node(parent->val*ratio, step, i, j);
                        parent->children.push_back(land_grid[i][j]);
                        land_mask[i][j] = true;
                        for(int x = -1; x < 1; x++) {
                            for(int y = -1; y <= 1; y++) {
                                if((x != 0 || y != 0) && i+x >= 0 && i+x < size && j+y >= 0 && j+y < size) {
                                    land_mask[i+x][j+y] = true;
                                }
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
                                    if((x != 0 || y != 0) && i+x >= 0 && i+x < size && j+y >= 0 && j+y < size) {

                                        neighborhood[inv_index] = land_grid[i+x][j+y];
                                        inv[inv_index] = p*neighborhood[inv_index]->val/(p*neighborhood[inv_index]->val+land_grid[i][j]->val*(1-p));
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
            }

            for(const auto& [i, j, val] : updates) {
                land_grid[i][j] = val;
                bool unmask = true;
                for(int x = -1; x < 1; x++) {
                    for(int y = -1; y <= 1; y++) {
                        if((x != 0 || y != 0) && i+x >= 0 && i+x < size && j+y >= 0 && j+y < size) {
                            bool unmask = true;
                            if(land_grid[i+x][j+x] != val) { unmask = false; }
                            for(int xx = -1; xx < 1; xx++) {
                                for(int yy = -1; yy <= 1; yy++) {
                                    if((xx != 0 || yy != 0) && i+x+xx >= 0 && i+x+xx < size && j+y+yy >= 0 && j+y+yy < size) {
                                        if(land_grid[i+x+xx][j+y+yy] != val) {
                                            unmask = false;
                                            break;
                                        }
                                    }
                                }
                                if(!unmask) break;
                            }
                            land_mask[i+x][j+y] = !unmask;
                        }
                    }
                }
            }
        }

        rep_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> rep_time = rep_end_time - rep_start_time;
        println("Rep %d run time = %fs", rep, rep_time.count());
        fflush(stdout);

        out_start_time = std::chrono::system_clock::now();
        FILE *fp;
        outfile = argv[3] + std::to_string(rep) + ".log";
        fp = fopen(outfile.c_str(), "w");
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size-1; j++) {
                fprintf(fp, "%f", land_grid[i][j]->val);
                fprintf(fp, ", ");
            }
            fprintf(fp, "%f\n", land_grid[i][size-1]->val);
        }
        fclose(fp);

        outfile = argv[3] + std::to_string(rep) + ".tree";
        fp = fopen(outfile.c_str(), "w");
        string output_tree = toString(speciation_root);
        fprintf(fp, "%s;", output_tree.c_str());
        println("Depth of tree = %d", get_depth(speciation_root));
        fclose(fp);
        delete_tree(speciation_root);

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

