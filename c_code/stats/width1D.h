// header file for width
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
#include "slope.h"
#include <time.h>
#include <chrono>
using namespace std;


// width, slope of width
tuple<array<double,10>,double> width1D(string fname) {
    int p2_start = 1;
    int p2_end = 10;
    int p2_range = p2_end-p2_start+1;
    // load landscape
    array<double,10> width_arr;
    fstream file (fname, ios::in);
    vector<double> landscape;
    string line, word; 
    int size = 0;
    if(file.is_open()){
        while(getline(file,line)){
            if (line == ""){
                continue;
            }
            stringstream str(line);
            while (getline(str,word,',')){
                landscape.push_back(log(stod(word)));
                size++;
            }
        }
    }
    else {
        cout<<"Could not open the file\n";
    }
    cout << size << endl;
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
            tempmean += landscape[i];
        }
        tempmean /= samp_size;
        for (int i=0; i<samp_size; ++i){
            tempvar += (landscape[i]-tempmean) * (landscape[i]-tempmean);
        }
        tempvar /= samp_size;
        double window_mean = tempmean;
        double window_var = tempvar;
        for (int i=samp_size; i<size; i++){
            double window_mean_old = window_mean;
            window_mean = window_mean_old + (landscape[i]-landscape[i-samp_size])/samp_size;
            window_var = window_var + pow(window_mean_old,2) - pow(window_mean,2) + (pow(landscape[i],2)-pow(landscape[i-samp_size],2))/samp_size;
            tempvar += window_var;
        }
        
        tempvar /= (size-samp_size+1);
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
    cout << "slope" << coeff << endl;
    return make_tuple(width_arr,coeff);
}










// simple width
double* width1D_simple(string fname) {
    int p2_start = 1;
    int p2_end = 20;
    int p2_range = p2_end-p2_start+1;
    // load landscape (landscape2 for conflict in landscape names)
    fstream file2 (fname, ios::in);
    vector<double> landscape2;
    string line2, word2; 
    int size = 0;
    if(file2.is_open()){
        while(getline(file2,line2)){
            if (line2 == ""){
                continue;
            }
            stringstream str(line2);
            while (getline(str,word2,',')){
                landscape2.push_back(log(stod(word2)));
                size++;
            }
        }
    }
    else {
        cout<<"Could not open the file\n";
    }
    // set<double> distinct_species;
    // for (int i=0; i<size; i++){
    //         distinct_species.insert(landscape2[i]);
    // }
    // cout << "diversity = " << distinct_species.size() << endl;
    // cout << size << endl;
    static double width_arr[20];
    vector<double> sample_list;
    for (int p2=p2_start; p2<=p2_end; ++p2){
        int samp_size = pow(2,p2);
        // sample_list.push_back(log(samp_size));
        // cout << "p2 = " << p2 << endl;
        // initialize window
        double tempwidth = 0;
        int nsamp = (int) size/samp_size;
        // cout << samp_size << endl;
        if (nsamp > 256 || samp_size==size){
        // if (true){
            for (int i=0; i<=nsamp-1; i++){
                double tempmean = 0;
                for (int ii=0; ii<samp_size; ii++){
                    tempmean += landscape2[i*samp_size+ii];
                }
                tempmean /= samp_size;
                for (int ii=0; ii<samp_size; ii++){
                    tempwidth += (landscape2[i*samp_size+ii]-tempmean)*(landscape2[i*samp_size+ii]-tempmean)/samp_size;
                }
            }
            tempwidth /= nsamp;
            width_arr[p2-p2_start] = tempwidth;
            // cout << "nsamp = " << nsamp << endl;
        }
        else {
            int nsamp_new = 256;
            int interval = (int) size/nsamp_new;
            for (int i=0; i<=nsamp_new-1; i++){
                double tempmean = 0;
                for (int ii=0; ii<samp_size; ii++){
                    tempmean += landscape2[(i*interval+ii)%size];
                }
                tempmean /= samp_size;
                for (int ii=0; ii<samp_size; ii++){
                    tempwidth += (landscape2[(i*interval+ii)%size]-tempmean)*(landscape2[(i*interval+ii)%size]-tempmean)/samp_size;
                }
            }
            tempwidth /= nsamp_new;
            width_arr[p2-p2_start] = tempwidth;
            // cout << "nsamp_new = " << nsamp_new << endl;
        }
        // cout << "samp_size = " << samp_size << "\twidth = " << tempwidth << endl;
    }
    return width_arr;
}


// simple width, parallel
double* width1D_simple_parallel(string fname) {
    std::chrono::time_point<std::chrono::system_clock> load_start_time, load_end_time, width_L_start_time, width_L_end_time; 
    int size = 1048576;
    double landscape[size];
    fstream file (fname, ios::in);
    load_start_time = std::chrono::system_clock::now();
    if (GA_Nodeid()==0){
        string line;
        getline(file,line);
        stringstream str(line);
        for (int i=0; i<size; i++){
            string word;
            getline(str, word, ',');
            landscape[i] = log(stod(word));
        }
    }
    MPI_Bcast(&landscape, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    load_end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> load_time = load_end_time - load_start_time;
    if (GA_Nodeid()==0){
        cout << "load time = " << load_time.count() << endl;
    }
    // cout << size << endl;
    // if (GA_Nodeid()==1){
    //     set<double> distinct_species;
    //     for (int i=0; i<size; i++){
    //             distinct_species.insert(landscape[i]);
    //     }
    //     cout << "diversity = " << distinct_species.size() << endl;
    // }
    // if (GA_Nodeid()==1){
    //     int boundary_count = 0;
    //     for (int i=0; i<size; i++){
    //         if (landscape[i] != landscape[(i-1+size)%size] || landscape[i] != landscape[(i+1)%size]){
    //             boundary_count += 1;
    //         }
    //     }
    //     cout << "actual landscape: number of cell on boundary = " << boundary_count << endl;
    // }
    int p2_start = 1;
    int p2_end = 20;
    int p2_range = p2_end-p2_start+1;
    static double width_arr[20];
    vector<double> sample_list;
    for (int p2=p2_start; p2<=p2_end; ++p2){
        int samp_size = pow(2,p2);
        double tempwidth = 0;
        int nsamp = (int) size/samp_size;
        width_L_start_time = std::chrono::system_clock::now();
        if (nsamp > 256 || samp_size==size){
        // if (true){
            for (int i=0; i<=nsamp-1; i++){
                double tempmean = 0;
                for (int ii=0; ii<samp_size; ii++){
                    tempmean += landscape[i*samp_size+ii];
                }
                tempmean /= samp_size;
                for (int ii=0; ii<samp_size; ii++){
                    tempwidth += (landscape[i*samp_size+ii]-tempmean)*(landscape[i*samp_size+ii]-tempmean)/samp_size;
                }
            }
            tempwidth /= nsamp;
            width_arr[p2-p2_start] = tempwidth;
            // cout << "nsamp = " << nsamp << endl;
        }
        else {
            int nsamp_new = 256;
            int interval = (int) size/nsamp_new;
            int nsamp_per_node = nsamp_new/GA_Nnodes();
            int i = GA_Nodeid();
            for (int samp_index=0; samp_index<nsamp_per_node; samp_index++){
                double tempmean = 0;
                for (int ii=0; ii<samp_size; ii++){
                    tempmean += landscape[((i*nsamp_per_node+samp_index)*interval+ii)%size];
                }
                tempmean /= samp_size;
                for (int ii=0; ii<samp_size; ii++){
                    tempwidth += pow(landscape[((i*nsamp_per_node+samp_index)*interval+ii)%size]-tempmean,2)/samp_size;
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, &tempwidth, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            tempwidth /= nsamp_new;
            width_arr[p2-p2_start] = tempwidth;
        }
        width_L_end_time = std::chrono::system_clock::now();
        std::chrono::duration<double> width_L_time = width_L_end_time - width_L_start_time;
        // if (GA_Nodeid()==0){
        //     cout << "samp_size = " << samp_size << "\t time = " << width_L_time.count() << "\twidth = " << tempwidth << endl;
        // }
    }
    return width_arr;
}