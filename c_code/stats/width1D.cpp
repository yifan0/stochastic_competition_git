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
#include "width1D.h"
using namespace std;


int main() {
    string fname = "1D_landscape_1024_5_40X/test_ga_1048576_1_128_rep0.csv";
    double *p = width1D_simple(fname);
}