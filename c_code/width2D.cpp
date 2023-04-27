// measure 2D slope of width for reps
#include <cstdio>
#include <random> // rand()
#include <string.h> // memcpy()
#include <queue> // queue
#include <iostream>
#include <fstream>
#include <omp.h>
#include <array>
#include <algorithm>
#include <vector>
#include <numeric>
#include "slope.h"
#include "width2D.h"
// #include "species_count.h"
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

int main() {
    ofstream myfile;
    vector<double> coeff_step;
    int nrep = 1;
    double coeff;
    double width_step[nrep][4];
    vector<int> species_count_step[nrep];
    vector<double> species_fitness_step[nrep];
    for (int rep=1; rep<=nrep; ++rep){
        cout << rep << endl;
        string fname = "test_ga_500_1_121_rep0.csv";
        tuple<array<double,4>,vector<int>,vector<double>,double> result = width2D(fname);
        // width
        array<double,4> width_arr = get<0>(result);
        // species count
        vector<int> species_count = get<1>(result);
        species_count_step[rep-1] = species_count;
        // species fitness
        vector<double> species_fitness = get<2>(result);
        species_fitness_step[rep-1] = species_fitness;
        // slope of width
        coeff = get<3>(result);
        coeff_step.push_back(coeff);
        for (int i=0; i<=3; i++){
            width_step[rep-1][i] = width_arr[i];
        }
    }
    myfile.open ("2D_width_500_4.csv");
    for (int i=0;i<=nrep-1;i++){
        for (int j=0; j<=3; j++){
            myfile << width_step[i][j] << ",";
        }
        myfile << endl;
    }
    myfile.close();
    myfile.open ("2D_species_count_500_4.csv");
    for (int i=0;i<=nrep-1;i++){
        vector<int> species_count = species_count_step[i];
        for (int count: species_count){
            myfile << count << ",";
        }
        myfile << endl;
    }
    myfile.close();
    myfile.open ("2D_species_fitness_500_4.csv");
    for (int i=0;i<=nrep-1;i++){
        vector<double> species_fitness = species_fitness_step[i];
        for (double fitness: species_fitness){
            myfile << fitness << ",";
        }
        myfile << endl;
    }
    myfile.close();    
    myfile.open ("2D_slope_of_width_500_8_4.csv");
    for (int i = 0; i <= nrep-1; i++){
        myfile << coeff_step[i] << ",";
    }
    myfile.close();
}