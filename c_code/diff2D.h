// 2D difference and diversity measure
// caveat: slow when there are many species (high specrate)
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

// difference, diversity
tuple<array<double,7>,array<double,7>> diff2D(string fname) {
    // load landscape
    array<double,7> diff_arr;
    array<double,7> div_arr;
    fstream file (fname, ios::in);
    vector<vector<double>> landscape;
    string line, word; 
    int size = 0;
    // int nrep = 0;
    if(file.is_open()){
        while(getline(file,line)){
            if (line == ""){
                // nrep += 1;
                // cout << landscape.size() << endl;
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
    vector<double> diff;
    vector<double> div;
    vector<double> sample_list;
    for (int p2=1; p2<=7; ++p2){
        int samp_size = pow(2,p2);
        sample_list.push_back(log(samp_size));
        cout << "p2 = " << p2 << endl;
        // initialize window
        double tempdiv = 0;
        set<double> distinct_element;
        for (int i=0; i<samp_size; ++i){
            for (int j=0; j<samp_size; ++j){
                distinct_element.insert(landscape[i][j]);
            }
        }
        tempdiv = distinct_element.size();
        int window_div = distinct_element.size();
        auto minmax = minmax_element(distinct_element.begin(),distinct_element.end());
        double tempmin = *minmax.first;
        double tempmax = *minmax.second;
        double window_min = tempmin;
        double window_max = tempmax;
        double tempdiff = tempmax - tempmin;
        double window_diff = tempdiff;
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
            set<double> toDrop;
            set<double> realDrop;
            set<double> realAdd;
            set<double> toAdd;
            // x even
            if (x%2==0){
                if (y==size-1){
                    for (int j=y; j<=y+samp_size-1; ++j){
                        toDrop.insert(landscape[x][j%size]);
                    }
                    for (set<double>::iterator ptr=toDrop.begin();ptr!=toDrop.end(); ++ptr){
                        realDrop.insert(*ptr);
                        bool found = false;
                        for (int i=x+1; i<=x+samp_size-1; ++i){
                            for (int j=y; j<=y+samp_size-1; ++j){
                                if (*ptr == landscape[i%size][j%size]){
                                    realDrop.erase(*ptr);
                                    found = true;
                                    break;
                                }
                            }
                            if (found){
                                break;
                            }
                        }
                    }
                    for (int j=y; j<=y+samp_size-1; ++j){
                        toAdd.insert(landscape[(x+samp_size)%size][j%size]);
                    }
                    for (set<double>::iterator ptr=toAdd.begin();ptr!=toAdd.end(); ++ptr){
                        realAdd.insert(*ptr);
                        bool found = false;
                        for (int i=x+samp_size-1; i>=x+1; --i){
                            for (int j=y; j<=y+samp_size-1; ++j){
                                if (*ptr == landscape[i%size][j%size]){
                                    realAdd.erase(*ptr);
                                    found = true;
                                    break;
                                }
                            }
                            if (found){
                                break;
                            }
                        }                            
                    }
                }
                else {
                    for (int i=x; i<=x+samp_size-1; ++i){
                        toDrop.insert(landscape[i%size][y]);
                    }
                    for (set<double>::iterator ptr=toDrop.begin();ptr!=toDrop.end(); ++ptr){
                        realDrop.insert(*ptr);
                        bool found = false;
                        for (int j=y+1; j<=y+samp_size-1; ++j){
                            for (int i=x; i<=x+samp_size-1; ++i){
                                if (*ptr == landscape[i%size][j%size]){
                                    realDrop.erase(*ptr);
                                    found = true;
                                    break;
                                }
                            }
                            if (found){
                                break;
                            }
                        }
                    }
                    for (int i=x; i<=x+samp_size-1; ++i){
                        toAdd.insert(landscape[i%size][(y+samp_size)%size]);
                    }
                    for (set<double>::iterator ptr=toAdd.begin();ptr!=toAdd.end(); ++ptr){
                        realAdd.insert(*ptr);
                        bool found = false;
                        for (int j=y+samp_size-1; j>=y+1; --j){
                            for (int i=x; i<=x+samp_size-1; ++i){
                                if (*ptr == landscape[i%size][j%size]){
                                    realAdd.erase(*ptr);
                                    found = true;
                                    break;
                                }
                            }
                            if (found){
                                break;
                            }
                        }                            
                    }
                }
            }
            // x odd
            else {
                if (y==0){
                    for (int j=y; j<=y+samp_size-1; ++j){
                        toDrop.insert(landscape[x][j]);
                    }
                    for (set<double>::iterator ptr=toDrop.begin();ptr!=toDrop.end(); ++ptr){
                        realDrop.insert(*ptr);
                        bool found = false;
                        for (int i=x+1; i<=x+samp_size-1; ++i){
                            for (int j=y; j<=y+samp_size-1; ++j){
                                if (*ptr == landscape[i%size][j%size]){
                                    realDrop.erase(*ptr);
                                    found = true;
                                    break;
                                }
                            }
                            if (found){
                                break;
                            }
                        }
                    }
                    for (int j=y; j<=y+samp_size-1; ++j){
                        toAdd.insert(landscape[(x+samp_size)%size][j]);
                    }
                    for (set<double>::iterator ptr=toAdd.begin();ptr!=toAdd.end(); ++ptr){
                        realAdd.insert(*ptr);
                        bool found = false;
                        for (int i=x+samp_size-1; i>=x+1; --i){
                            for (int j=y; j<=y+samp_size-1; ++j){
                                if (*ptr == landscape[i%size][j]){
                                    realAdd.erase(*ptr);
                                    found = true;
                                    break;
                                }
                            }
                            if (found){
                                break;
                            }
                        }
                    }
                }
                else {
                    for (int i=x; i<=x+samp_size-1; ++i){
                        toDrop.insert(landscape[i%size][(y+samp_size-1)%size]);
                    }
                    for (set<double>::iterator ptr=toDrop.begin();ptr!=toDrop.end(); ++ptr){
                        realDrop.insert(*ptr);
                        bool found = false;
                        for (int j=y+samp_size-2; j>=y; j--){
                            for (int i=x; i<=x+samp_size-1; ++i){
                                if (*ptr == landscape[i%size][j%size]){
                                    realDrop.erase(*ptr);
                                    found = true;
                                    break;
                                }
                            }
                            if (found){
                                break;
                            }
                        }
                    }
                    for (int i=x; i<=x+samp_size-1; ++i){
                        toAdd.insert(landscape[i%size][y-1]);
                    }
                    for (set<double>::iterator ptr=toAdd.begin();ptr!=toAdd.end(); ++ptr){
                        realAdd.insert(*ptr);
                        bool found = false;
                        for (int j=y; j<=y+samp_size-2; ++j){
                            for (int i=x; i<=x+samp_size-1; ++i){
                                if (*ptr == landscape[i%size][j%size]){
                                    realAdd.erase(*ptr);
                                    found = true;
                                    break;
                                }
                            }
                            if (found){
                                break;
                            }
                        }
                    }
                }
            }
            for (set<double>::iterator ptr=realAdd.begin();ptr!=realAdd.end(); ++ptr){
                if (*ptr > tempmax){
                    tempmax = *ptr;
                }
                if (*ptr < tempmin){
                    tempmin = *ptr;
                }
            }
            // order
            if (realDrop.count(tempmin) || realDrop.count(tempmax)){
                set<double> distinct_sample;
                if (x%2==0){
                    if (y==size-1){
                        for (int i=x+1; i<=x+samp_size; i++){
                            for (int j=y; j<=y+samp_size-1; j++){
                                distinct_sample.insert(landscape[i%size][j%size]);
                            }
                        }
                    }
                    else {
                        for (int i=x; i<=x+samp_size-1; i++){
                            for (int j=y+1; j<=y+samp_size; j++){
                                distinct_sample.insert(landscape[i%size][j%size]);
                            }
                        }
                    }
                }
                else {
                    if (y==0){
                        for (int i=x+1; i<=x+samp_size; i++){
                            for (int j=y; j<=y+samp_size-1; j++){
                                distinct_sample.insert(landscape[i%size][j%size]);
                            }
                        }
                    }
                    else {
                        for (int i=x; i<=x+samp_size-1; i++){
                            for (int j=y-1; j<=y+samp_size-2; j++){
                                distinct_sample.insert(landscape[i%size][j%size]);
                            }
                        }
                    }
                }
                auto minmax = minmax_element(distinct_sample.begin(),distinct_sample.end());
                tempmin = *minmax.first;
                tempmax = *minmax.second;
            }
            window_div += (realAdd.size() - realDrop.size());
            tempdiv += window_div;
            window_diff = tempmax - tempmin;
            tempdiff += window_diff;
        }
        tempdiv /= (size*size);
        // cout << "tempdiv" << tempdiv << endl;
        // cout << log(tempdiv) << endl;
        tempdiff /= (size*size);
        // cout << "tempdiff" << tempdiff << endl;
        // cout << log(tempdiff) << endl;
        // if (tempvar == 0 || log(tempvar) < -10){
        //     cout << endl;
        //     return make_tuple(width_arr,0.0);
        // }
        // for calculation
        // diff.push_back(log(tempdiff));
        // div.push_back(log(tempdiv));
        // for record
        diff_arr[p2-1] = tempdiff;
        div_arr[p2-1] = tempdiv;
    }
    cout << endl;
    // double coeff_diff = slope(sample_list,diff);
    // double coeff_div = slope(sample_list,div);
    return make_tuple(diff_arr,div_arr);
}