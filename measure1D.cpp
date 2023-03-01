// move window -- diversity & difference & width
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <typeinfo>
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

int main() {
    string fname = "landscape1D512_10.csv";
    fstream file (fname, ios::in);
    vector<double> landscape;
    string line, word; 
    int size = 0;
    if(file.is_open()){
        while(getline(file,word,',')){
            landscape.push_back(stod(word));
            size++;
        }
    }
    else {
        cout<<"Could not open the file\n";
    }
    // vector<double> landscape;
    // int size = 20;
    // for (int i=0;i<size;i++){
    //     landscape.push_back(round(i/5));
    //     // landscape.push_back(i%10);
    // } 
    // cout << "landscape:" << endl;
    // for (int i=0; i<landscape.size();i++){
    //     if (i==10){
    //         cout << "\n";
    //     }
    //     cout << "\t" << landscape[i];
    // }
    // window of size 10
    vector<double> diversity;
    vector<double> difference;
    vector<double> width;
    for (int p2=1; p2<18; p2++){
        int samp_size = pow(2,p2);
        cout << p2 << endl;
        // int samp_size = 5;
        // initialize window
        vector<double> samp;
        for (int i=0; i<samp_size; i++){
            samp.push_back(landscape[i]);
        }
        // calculation on the initial window
        double tempdiv = 0;
        sort(samp.begin(),samp.end());
        auto it = unique(samp.begin(),samp.end());
        samp.resize(distance(samp.begin(), it));
        int window_div = samp.size();
        tempdiv += window_div;
        auto result = minmax_element(samp.begin(),samp.end());
        double tempmin = *result.first;
        double tempmax = *result.second;
        double window_min = tempmin;
        double window_max = tempmax;
        double tempdiff = 0;
        double window_diff = tempmax - tempmin;
        tempdiff += window_diff;
        double tempwidth = 0;
        double tempmean = 0;
        double tempvar = 0;
        for (int i=0; i<samp_size; i++){
            tempmean += landscape.at(i);
        }
        tempmean /= samp_size;
        for (int i=0; i<samp_size; i++){
            tempvar += (landscape.at(i)-tempmean) * (landscape.at(i)-tempmean);
        }
        tempvar /= samp_size;
        double window_mean = tempmean;
        double window_var = tempvar;
        // cout << window_mean << endl;
        // scroll window
        // only look at the value, avoid push_back
        // sample from [i-samp_size+1] to [i], new element [i]
        // start moving, i starts from samp_size
        // old sample [i-samp_size],...,[i-1]
        // new sample [i-samp_size+1],...,[i]
        // sample overlap [i-samp_size+1],...[i-1]
        for (int i=samp_size; i<size; i++){
            // cout << i-samp_size+1 << "," << i << endl;
            // iterate over sample overlap, search for [i-samp_size] and [i]
            bool drop_change = true;
            bool add_change = true;
            for (int j=i-samp_size+1; j<i; j++){
                // [i-samp_size] exists in the overlap, ok to drop
                if (landscape.at(j) == landscape.at(i-samp_size)){
                    drop_change = false;
                    break;
                }
            }
            // for (int j=i-samp_size+1; j<i; j++){
            for (int j=i-1; j>i-samp_size; j--){
                // [i] exists in the overlap, no need to add
                if (landscape.at(j) == landscape.at(i)){
                    add_change = false;
                    break;
                }
            }
            window_div += (add_change - drop_change);
            tempdiv += window_div;
            // cout << "\tdrop" << drop_change;
            // cout << "\tadd" << add_change;
            // difference
            if (drop_change){
                if (landscape.at(i-samp_size) == tempmin){
                    tempmin = *min_element(landscape.begin()+i-samp_size+1,landscape.begin()+i);
                } 
                if (landscape.at(i-samp_size) == tempmax){
                    tempmax = *max_element(landscape.begin()+i-samp_size+1,landscape.begin()+i);
                } 
            } 
            if (add_change){
                if (landscape.at(i) < tempmin){
                    tempmin = landscape.at(i);
                } 
                if (landscape.at(i) > tempmax){
                    tempmax = landscape.at(i);
                } 
            }
            window_diff = tempmax - tempmin;
            tempdiff += window_diff;
            // cout << "\twindiff" << window_diff << endl;
            // width
            // recalculate window_var 
            double window_mean_old = window_mean;
            window_mean = window_mean_old + (landscape.at(i) - landscape.at(i-samp_size))/samp_size;
            window_var = window_var + pow(window_mean_old,2) - pow(window_mean, 2) + (pow(landscape.at(i),2)-pow(landscape.at(i-samp_size),2))/samp_size;
            tempvar += window_var;
            // cout << "\tmean" << window_mean;
            // cout << "\tvar" << window_var << endl;
        }
        // cout << "\ntempdiff" << tempdiff << endl;
        tempdiv /= (size-samp_size+1);
        tempdiff /= (size-samp_size+1);
        tempvar /= (size-samp_size+1);
        // cout << tempvar << endl;
        diversity.push_back(tempdiv);
        difference.push_back(tempdiff);
        width.push_back(tempvar);
    }
    ofstream myfile;
    myfile.open ("test_div.csv");
    for (int i = 0; i < diversity.size()-1; i++){
        myfile << diversity[i] << ",";
    }
    myfile << diversity[diversity.size()-1];
    myfile.close();
    myfile.open ("test_diff.csv");
    for (int i = 0; i < difference.size()-1; i++){
        myfile << difference[i] << ",";
    }
    myfile << difference[difference.size()-1];
    myfile.close();
    myfile.open ("test_width.csv");
    for (int i = 0; i < width.size()-1; i++){
        myfile << width[i] << ",";
    }
    myfile << width[width.size()-1];
    myfile.close();
}
