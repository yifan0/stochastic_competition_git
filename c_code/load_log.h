// convert GA output log file to csv
#include <cstdio>
#include <random> // rand()
#include <string.h> // memcpy()
#include <queue> // queue
#include <iostream>
#include <fstream>
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
#include <set>
#include "width1D.h"
#include <typeinfo>
using namespace std;


void log_to_csv(string fname) {
    fstream file (fname, ios::in);
    vector<double> landscape;
    set<double> distinct_species;
    string line, word; 
    int size = 0;
    // int nrep = 0;
    int counter = 0;
    int count_line = 0;
    if(file.is_open()){
        while(getline(file,line)){
            count_line += 1;
            if (line == ""){
                // skip empty line
                continue;
            }
            stringstream str(line); 
            vector<double> row;
            while (getline(str, word, ' ')){
                if (word != "" && count_line > 2){
                    counter += 1;
                    // cout << "counter = " << counter << endl;
                    if (counter %2 == 0){
                        // cout << word << endl;
                        landscape.push_back(stod(word));
                        distinct_species.insert(stod(word));
                        size += 1;
                    }
                }
            }
        }
    }
    else {
        cout<<"log_to_csv could not open the file\n";
    }
    // println("size = %d", size);
    // println("diversity = %d", distinct_species.size());
    string fname_new =  fname.substr(0,fname.length()-4)+".csv";
    ofstream myfile;
    double mean = 0;
    for (int i=0; i<size-1; i++){
            mean += landscape[i]/size;
    }
    // println("log-csv mean = %d", mean);
    // cout << "log-csv mean = " << mean << endl;
    myfile.open (fname_new.c_str());
    for (int i=0; i<size-1; i++){
            myfile << landscape[i] << ",";
    }
    myfile << landscape[size-1];
    myfile.close();
}