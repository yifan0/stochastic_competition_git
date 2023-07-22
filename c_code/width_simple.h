// measure 2D slope of width for reps, samples taken side by side
// functions width_simple and width_simple_parallel
#include <cstdio>
#include <random> // rand()
#include <string.h> // memcpy()
#include <queue> // queue
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <array>
#include <algorithm>
#include <vector>
#include <numeric>
#include <time.h>
#include <chrono>
using namespace std;

// 2^1, 2^1.5, 2^2, ..., 2^11
double* width_simple(string fname) {
    // std::chrono::time_point<std::chrono::system_clock> width_L_start_time, width_L_end_time;
    // string fname = "2D_landscape_2048_7/test_ga_2048_1_121_rep0_checkpoint14.csv";
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
    
    int start_pt = 2;
    int end_pt = ((int) log2(size))*2;
    // cout << end_pt << endl;
    static double width_arr[23];
    vector<double> sample_list;
    for (int p2=start_pt; p2<=end_pt; ++p2){
        int samp_size = (int) pow(sqrt(2),p2);
        // cout << "samp_size = " << samp_size << endl;
        // sample_list.push_back(log(samp_size));
        // cout << "p2 = " << p2 << endl;
        // width_L_start_time = std::chrono::system_clock::now();
        int nsamp = size/samp_size;
        // if number of samples is > 256 or sample is the whole landscape
        if (nsamp > 16 || samp_size == size){
            double tempdiv = 0;
            double tempdiff = 0;
            double tempwidth = 0;
            for (int i=0; i<=nsamp-1; i++){
                for (int j=0; j<=nsamp-1; j++){
                    double tempmean = 0;
                    for (int ii=0; ii<samp_size; ++ii){
                        for (int jj=0; jj<samp_size; ++jj){
                            tempmean += landscape[i*samp_size+ii][j*samp_size+jj];
                        }
                    }
                    tempmean /= (samp_size*samp_size);
                    for (int ii=0; ii<samp_size; ++ii){
                        for (int jj=0; jj<samp_size; ++jj){
                            tempwidth += pow(landscape[i*samp_size+ii][j*samp_size+jj]-tempmean,2)/(samp_size*samp_size);
                        }
                    }
                }
            }
            tempwidth /= (nsamp*nsamp);
            width_arr[p2-start_pt] = tempwidth;
        }
        // if number of samples < 256, manually make it than 256
        else {
            double tempdiv = 0;
            double tempdiff = 0;
            double tempwidth = 0;
            int nsamp_new = 16;
            // interval/gap between samples
            int interval = (int) size/nsamp_new;
            for (int i=0; i<nsamp_new; i++){
                for (int j=0; j<nsamp_new; j++){
                    double tempmean = 0;
                    for (int ii=0; ii<samp_size; ++ii){
                        for (int jj=0; jj<samp_size; ++jj){
                            tempmean += landscape[(i*interval+ii)%size][(j*interval+jj)%size];
                        }
                    }
                    tempmean /= (samp_size*samp_size);
                    for (int ii=0; ii<samp_size; ++ii){
                        for (int jj=0; jj<samp_size; ++jj){
                            tempwidth += pow(landscape[(i*interval+ii)%size][(j*interval+jj)%size]-tempmean,2)/(samp_size*samp_size);
                        }
                    }
                }
            }
            tempwidth /= (nsamp_new*nsamp_new);
            width_arr[p2-start_pt] = tempwidth;
        }
        // width_L_end_time = std::chrono::system_clock::now();
        // std::chrono::duration<double> width_L_time = width_L_end_time - width_L_start_time;
        // cout << "p2 = " << p2 << "measure time = " << width_L_time.count() << endl;
    }
    return width_arr;
}

double* width_simple_parallel(string fname) {
    // string fname = "2D_landscape_2048_7/test_ga_2048_1_121_rep0_checkpoint14.csv";
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
    
    int start_pt = 2;
    int end_pt = ((int) log2(size))*2;
    // cout << end_pt << endl;
    static double width_arr[23];
    vector<double> sample_list;
    for (int p2=start_pt; p2<=end_pt; ++p2){
        int samp_size = (int) pow(sqrt(2),p2);
        // cout << "samp_size = " << samp_size << endl;
        // sample_list.push_back(log(samp_size));
        // cout << "p2 = " << p2 << endl;
        int nsamp = size/samp_size;
        // if number of samples is > 256 or sample is the whole landscape
        if (nsamp > 16 || samp_size == size){
            double tempdiv = 0;
            double tempdiff = 0;
            double tempwidth = 0;
            for (int i=0; i<=nsamp-1; i++){
                for (int j=0; j<=nsamp-1; j++){
                    double tempmean = 0;
                    for (int ii=0; ii<samp_size; ++ii){
                        for (int jj=0; jj<samp_size; ++jj){
                            tempmean += landscape[i*samp_size+ii][j*samp_size+jj];
                        }
                    }
                    tempmean /= (samp_size*samp_size);
                    for (int ii=0; ii<samp_size; ++ii){
                        for (int jj=0; jj<samp_size; ++jj){
                            tempwidth += pow(landscape[i*samp_size+ii][j*samp_size+jj]-tempmean,2)/(samp_size*samp_size);
                        }
                    }
                }
            }
            tempwidth /= (nsamp*nsamp);
            width_arr[p2-start_pt] = tempwidth;
        }
        // if number of samples < 256, manually make it than 256
        else {
            double tempdiv = 0;
            double tempdiff = 0;
            double tempwidth = 0;
            int nsamp_new = 16;
            // interval/gap between samples
            int interval = (int) size/nsamp_new;
            // assign local area based on node id
            int i = GA_Nodeid()/16;
            int j = GA_Nodeid()%16;
            double tempmean = 0;
            for (int ii=0; ii<samp_size; ++ii){
                for (int jj=0; jj<samp_size; ++jj){
                    tempmean += landscape[(i*interval+ii)%size][(j*interval+jj)%size];
                }
            }
            tempmean /= (samp_size*samp_size);
            for (int ii=0; ii<samp_size; ++ii){
                for (int jj=0; jj<samp_size; ++jj){
                    tempwidth += pow(landscape[(i*interval+ii)%size][(j*interval+jj)%size]-tempmean,2)/(samp_size*samp_size);
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, &tempwidth, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            tempwidth /= GA_Nnodes();
            width_arr[p2-start_pt] = tempwidth;
        }

    }
    return width_arr;
}