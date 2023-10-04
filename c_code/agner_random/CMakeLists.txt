cmake_minimum_required (VERSION 3.18.0 FATAL_ERROR)
project (Agner VERSION 1.0.0 LANGUAGES CXX)

add_library(Agner randomc.h  sfmt.h  stocc.h  mother.cpp    stoc3.cpp  wnchyppr.cpp    rancombi.cpp
   sfmt.cpp
fnchyppr.cpp  stoc1.cpp mersenne.cpp  stoc2.cpp  userintf.cpp)
