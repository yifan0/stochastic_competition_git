#include <cstdio>
#include <random> // rand()
#include <string.h> // memcpy()
#include <queue> // queue
#include <chrono>
#include <set>
#include <map>
#include <numeric>
#include <algorithm>
#include "sfmt.h"
using namespace std;

#ifndef MULTIFILE_PROJECT
// If compiled as a single file then include these cpp files, 
// If compiled as a project then compile and link in these cpp files.
    // Include code for the chosen random number generator:
    #include "sfmt.cpp"
    // define system specific user interface:
    #include "userintf.cpp"
#endif

#define println(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

#define GA_MAX(a,b) (((a)>(b)) ? (a) : (b))
#define GA_MIN(a,b) (((a)<(b)) ? (a) : (b))

int main(int argc, char* argv[]) {
    std::chrono::time_point<std::chrono::system_clock> start_time, end_time;
    double p = 0.1;
    size_t size = 500;
    size_t n = size*size;
    size_t steps = (size*1.0/p);
    vector<int> vec(n);

    println("Size = %d, steps = %d", size, steps);

    // Random Number Generator
    CRandomSFMT1 RanGen((int)time(NULL));

    start_time = std::chrono::system_clock::now();

    for(int t = 0; t < steps; t++) {
        for(int i = 0; i < n; i++) {
            vec[i] = RanGen.Random();
        }
    }

    end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;
    println("RNG time = %f s", total_time.count());

    start_time = std::chrono::system_clock::now();
    size_t row, col;
    double max = 0;
    for(int t = 0; t < steps; t++) {
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                max = vec[i*size+j];
                for(int x = -1; x <= 1; x++) {
                    for(int y = -1; y <= 1; y++) {
                        if((x != 0 || y != 0)) {
                            row = i+x;
                            col = j+y;
                            if(i+x < 0 || i+x >= size)
                                row = (row+size)%size;
                            if(j+y < 0 || j+x >= size)
                                col = (col+size)%size;
                            max = GA_MAX(max, vec[row*size+col]);
                        }
                    }
                }
                vec[i*size+j] = max;
            }
        }
    }
    end_time = std::chrono::system_clock::now();
    total_time = end_time - start_time;
    println("Visiting neighbors with wrap-around time = %f s", total_time.count());

    start_time = std::chrono::system_clock::now();
    for(int t = 0; t < steps; t++) {
        for(int i = 0; i < n; i++) {
                max = vec[i];
                for(int x = -4; x <= 4; x++) {
                    if(x != 0 && i+x >= 0 && i+x < size) {
                        max = GA_MAX(max, vec[i+x]);
                    }
                }
                vec[i] = max;
        }
    }
    end_time = std::chrono::system_clock::now();
    total_time = end_time - start_time;
    println("Visiting neighbors single row time = %f s", total_time.count());

    return 0;
}

