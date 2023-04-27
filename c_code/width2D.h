// 2D width measure
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <tuple>
#include <math.h>
#include <string>
#include <sstream>
#include <typeinfo>
#include <set>
#include "species_count.h"
// #include "slope.h"
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

// double slope(const std::vector<double>& x, const std::vector<double>& y) {
//     const auto n    = x.size();
//     const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
//     const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
//     const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
//     const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
//     const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
//     return a;
// }
// width, species count, species fitness, slope of width
tuple<array<double,4>,vector<int>,vector<double>,double> width2D(string fname) {
    int p2_start = 4;
    int p2_end = 7;
    int p2_range = p2_end-p2_start+1;
    // load landscape
    array<double,4> width_arr;
    width_arr[0]=0.0;
    width_arr[1]=0.0;
    width_arr[2]=0.0;
    width_arr[3]=0.0;
    fstream file (fname, ios::in);
    vector<vector<double>> landscape;
    string line, word; 
    int size = 0;
    // int nrep = 0;
    if(file.is_open()){
        while(getline(file,line)){
            if (line == ""){
                // nrep += 1;
                cout << landscape.size() << endl;
                // skip empty line
                continue;
            }
            stringstream str(line); 
            vector<double> row;
            while (getline(str, word, ',')){
                row.push_back(log(stod(word)));
            }
            landscape.push_back(row);
        }
    }
    else {
        cout<<"Could not open the file\n";
    }
    size = landscape[0].size();
    tuple<int, vector<int>, vector<double>> result = species_count(size, landscape);
    int div = get<0>(result);
    vector<int> species_count = get<1>(result);
    vector<double> species_fitness = get<2>(result);
    if (div == 1){
        return make_tuple(width_arr,species_count, species_fitness, 0.0);
    }
    vector<double> width;
    vector<double> sample_list;
    for (int p2=p2_start; p2<=p2_end; ++p2){
        int samp_size = pow(2,p2);
        sample_list.push_back(log(samp_size));
        // cout << "p2 = " << p2 << endl;
        // initialize window
        double tempwidth = 0;
        double tempmean = 0;
        double tempvar = 0;
        for (int i=0; i<samp_size; ++i){
            for (int j=0; j<samp_size; ++j){
                tempmean += landscape[i][j];
            }
        }
        tempmean /= (samp_size*samp_size);
        for (int i=0; i<samp_size; ++i){
            for (int j=0; j<samp_size; ++j){
                tempvar += (landscape[i][j]-tempmean) * (landscape[i][j]-tempmean);
            }
        }
        tempvar /= (samp_size*samp_size);
        double window_mean = tempmean;
        double window_var = tempvar;
        // create a list for the sample window movement
        vector<tuple<int,int>> traj;
        for (int i=0; i<=size-1; ++i){
            if (i%2==0){
                for (int j=0; j<=size-1; ++j){
                    traj.push_back(tuple<int,int>{i,j});
                }
            }
            else {
                for (int j=size-1;j>=0;j--){
                    traj.push_back(tuple<int,int>{i,j});
                }
            }
        }
        // drop the last traj point
        traj.pop_back();
        // sample window moves along the trajectory traj
        // at traj[i], sample is [get<0>(traj[i]):get<0>(traj[i])+samp_size-1, get<1>(traj[i]):get<1>(traj[i])+samp_size-1]
        // overlapping region
        // case 1: get<0> even
        // special case 1: get<0> even && get<1> == size-samp_size
        // case 2: get<0> odd
        // special case 2: get<0> odd && get<1> == 0
        // position of indicator of previous sample
        for (auto [x,y]:traj){
            double window_mean_old = window_mean;
            double addToMean = 0;
            double addToVar = 0;
            // x even
            if (x%2==0){
                if (y==size-1){
                    for (int j=y; j<=y+samp_size-1; ++j){
                        addToMean -= landscape[x][j%size];
                        addToVar -= pow(landscape[x][j%size],2);
                    }
                    for (int j=y; j<=y+samp_size-1; ++j){
                        addToMean += landscape[(x+samp_size)%size][j%size];
                        addToVar += pow(landscape[(x+samp_size)%size][j%size],2);
                    }
                }
                else {
                    for (int i=x; i<=x+samp_size-1; ++i){
                        addToMean -= landscape[i%size][y];
                        addToVar -= pow(landscape[i%size][y],2);
                    }
                    for (int i=x; i<=x+samp_size-1; ++i){
                        addToMean += landscape[i%size][(y+samp_size)%size];
                        addToVar += pow(landscape[i%size][(y+samp_size)%size],2);
                    }
                }
            }
            // x odd
            else {
                if (y==0){
                    for (int j=y; j<=y+samp_size-1; ++j){
                        addToMean -= landscape[x][j];
                        addToVar -= pow(landscape[x][j],2);
                    }
                    for (int j=y; j<=y+samp_size-1; ++j){
                        addToMean += landscape[(x+samp_size)%size][j];
                        addToVar += pow(landscape[(x+samp_size)%size][j],2);
                    }
                }
                else {
                    for (int i=x; i<=x+samp_size-1; ++i){
                        addToMean -= landscape[i%size][(y+samp_size-1)%size];
                        addToVar -= pow(landscape[i%size][(y+samp_size-1)%size],2);
                    }
                    for (int i=x; i<=x+samp_size-1; ++i){
                        addToMean += landscape[i%size][y-1];
                        addToVar += pow(landscape[i%size][y-1],2);
                    }
                }
            }
            window_mean = window_mean_old + addToMean/(samp_size*samp_size);
            window_var = window_var + pow(window_mean_old,2) - pow(window_mean,2) + addToVar/(samp_size*samp_size);
            tempvar += window_var;
        }
        tempvar /= (size*size);
        // cout << "var" << tempvar << endl;
        // cout << "log var" << log(tempvar) << endl;
        cout << log(tempvar) << "\t";
        // if (tempvar == 0 || log(tempvar) < -10){
        //     cout << endl;
        //     return make_tuple(width_arr,0.0);
        // }
        width.push_back(log(tempvar));
        width_arr[p2-p2_start] = tempvar;
    }
    cout << endl;
    double coeff = slope(sample_list,width);
    // cout << "slope" << coeff << endl;
    return make_tuple(width_arr,species_count,species_fitness,coeff);
}