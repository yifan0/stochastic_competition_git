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
    int size = 5;
    int nrep = 10;
    // double landrept[size][size];
    vector<vector<double>> landrept;
    for (int i=0; i<size; i++){
        vector<double> v1;
        for (int j=0; j<size; j++){
            // landrept[i][j] = i+j;
            // landrept.push_back((i+j));
            v1.push_back(i+j);
        }
        landrept.push_back(v1);
    }
    // string fname = "/home/yifan12/stochastic_competition_git/c_code/cpptest_1k_grid_result.txt";
    // vector<vector<double>> content;
    // int size = 0;
    // string line, word; 
    // fstream file (fname, ios::in);
    // if(file.is_open()){
    //     while(getline(file, line)){
    //     vector<double> row;
    //     stringstream str(line);
    //     while(getline(str, word, ','))
    //         row.push_back(stod(word));
    //         content.push_back(row);
    //         size++;
    //     }
    // }
    // else {
    //     cout<<"Could not open the file\n";
    // }
    // cout << size;
    double diversity[size][size][size];
    double width[size][size][size];
    double difference[size][size][size];
    for (int side=0; side < size; side++){
        for (int i=0; i<size-side; i++){
            for (int j=0; j<size-side; j++){
                vector<double> samp;
                for (int samp1=0; samp1<side+1; samp1++){
                    for (int samp2=0; samp2<side+1; samp2++){
                        samp.push_back(landrept[i+samp1][j+samp2]);
                    }
                }
                double mean = 0;
                for (int n=0; n<(side+1)*(side+1); n++){
                    mean += samp[n]/((side+1)*(side+1));
                }
                double var = 0;
                for (int n=0; n<((side+1)*(side+1)); n++){
                    var += (samp[n]-mean) * (samp[n]-mean);
                }
                width[side][i][j] = var/((side+1)*(side+1));
                // if (side == 2){
                //     cout<<"Printing samp"<<endl;
                //     for (double c : samp){
                //         cout<<c<<"\n";
                //     }
                // }
                // 
                // could use set.insert()
                sort(samp.begin(),samp.end());
                auto it = unique(samp.begin(),samp.end());
                samp.resize(distance(samp.begin(), it));
                diversity[side][i][j] = samp.size();
                // if (side == 2){
                //     cout<<"Printing unique"<<endl;
                //     for (double c : samp){
                //         cout<<c<<"\n";
                //     }
                // }
                // if (side == 2){
                //     cout<<"Printing diversity"<<endl;
                //     cout<<diversity[side][i][j]<<endl;
                // }
                // difference[side][i][j] = max(samp) - min(samp);
                auto result = minmax_element(samp.begin(),samp.end());
                difference[side][i][j] = *result.second - *result.first;
                // if (side==2){
                // auto result = minmax_element(samp.begin(),samp.end());
                // cout<<*result.second - *result.first<<endl;
                // }
            }
        }
    }
    cout<<"Printing landrept:\n";
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{
			cout<<"\t"<<landrept[i][j];
		}
		cout<<endl;
	}
    cout<<"Printing diversity:\n";
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{
			cout<<"\t"<<diversity[2][i][j];
		}
		cout<<endl;
	}
    cout<<"Printing width:\n";
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{
			cout<<"\t"<<width[2][i][j];
		}
		cout<<endl;
	}
    cout<<"Printing difference:\n";
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{
			cout<<"\t"<<difference[2][i][j];
		}
		cout<<endl;
	}
}