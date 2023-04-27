// measure 2D difference and diversity  for reps
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
#include "diff2D.h"
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }


int main() {
    ofstream myfile;
    int nrep = 1;
    double diff_step[nrep][7];
    double div_step[nrep][7];
    for (int rep=1; rep<=nrep; ++rep){
        cout << rep << endl;
        string fname = "test_ga_500_1_121_rep0.csv";
        tuple<array<double,7>,array<double,7>> result = diff2D(fname);
        // diff
        array<double,7> diff_arr = get<0>(result);
        // div
        array<double,7> div_arr = get<1>(result);
        for (int i=0; i<=6; i++){
            diff_step[rep-1][i] = diff_arr[i];
            div_step[rep-1][i] = div_arr[i];
        }
    }
    myfile.open ("2D_diff_500_4.csv");
    for (int i=0;i<=nrep-1;i++){
        for (int j=0; j<=6; j++){
            myfile << diff_step[i][j] << ",";
        }
        myfile << endl;
    }
    myfile.close();
    myfile.open ("2D_div_500_4.csv");
    for (int i=0;i<=nrep-1;i++){
        for (int j=0; j<=6; j++){
            myfile << div_step[i][j] << ",";
        }
        myfile << endl;
    }
    myfile.close();    
}