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
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

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










// width, slope of width
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
        if (nsamp > 128 || samp_size==size){
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
            int nsamp_new = 128;
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
        // cout << tempwidth << ",";
    }
    return width_arr;
}