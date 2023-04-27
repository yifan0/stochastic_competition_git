// species count
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
using namespace std;

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }
// div, species count, species fitness
tuple<int, vector<int>, vector<double>> species_count(int size, vector<vector<double>> land_grid) {
    set<double> distinct_element;
    for (int i=0; i<size; ++i){
        for (int j=0; j<size; ++j){
            distinct_element.insert(land_grid[i][j]);
        }
    }
    int div = distinct_element.size();
    vector<int> species_count;
    vector<double> species_fitness;
    cout << "div = " << div << endl;
    // cout << "species count" << endl;
    for (set<double>::iterator ptr = distinct_element.begin();ptr!=distinct_element.end();++ptr){
        // cout << *ptr << ":\t";
        int nCount = 0;
        for (int i=0;i<size;++i){
            nCount += count(land_grid[i].begin(),land_grid[i].end(),*ptr);
        }
        // cout << nCount << endl;
        species_fitness.push_back(*ptr);
        species_count.push_back(nCount);
    }
    // cout << endl;
    return make_tuple(div,species_count,species_fitness);
}