#ifndef PRINT_STATS_H
#define PRINT_STATS_H

#include <vector>
#include <string>
#include <fstream>
using namespace std;

int print_double_stat_to_csv(vector<double> stat, std::string file_name) {
    ofstream myfile;
    myfile.open(file_name);
    for (int i = 0; i < stat.size()-1; i++) {
        myfile << stat[i] << ",";
    }
    myfile << stat[stat.size()-1];
    myfile.close();
    return 0;
}

int print_double_arr_to_csv(double* stat, int len, std::string file_name) {
    ofstream myfile;
    myfile.open(file_name);
    for (int i = 0; i < len-1; i++) {
        myfile << stat[i] << ",";
    }
    myfile << stat[len-1];
    myfile.close();
    return 0;
}

int print_double_to_csv(double val, std::string file_name) {
    ofstream myfile;
    myfile.open(file_name);
    myfile << val;
    myfile.close();
    return 0;
}

int print_int_stat_to_csv(vector<int> stat, std::string file_name) {
    ofstream myfile;
    myfile.open(file_name);
    for (int i = 0; i < stat.size()-1; i++) {
        myfile << stat[i] << ",";
    }
    myfile << stat[stat.size()-1];
    myfile.close();
    return 0;
}

#endif // PRINT_STATS_H