cmake_minimum_required (VERSION 3.18.0 FATAL_ERROR)
project (stochastic_competition VERSION 1.0.0 LANGUAGES CXX)

find_package (MPI REQUIRED)
find_package(OpenMP)
#find_package(UPCXX REQUIRED)

set(ENABLE_FORTRAN OFF)

include(FetchContent)
FetchContent_Declare(
    globalarrays
    GIT_REPOSITORY https://github.com/kodingkoning/ga.git
)
FetchContent_Declare(
    mxx
    GIT_REPOSITORY https://github.com/patflick/mxx.git
)
FetchContent_MakeAvailable(globalarrays)
FetchContent_MakeAvailable(mxx)

include_directories(agner_random include)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++17")

set(targets
    sim
    omp_sim
    ga_sim
    ga_sim_no_tree
    ga_sim1D
    ga_sim1D_hydro
    benchmark
    tree_sim
)

foreach(TARGET ${targets})
    add_executable(${TARGET} "${TARGET}.cpp")
    target_compile_features(${TARGET} PUBLIC cxx_std_17)
endforeach(TARGET)

target_include_directories(ga_sim PUBLIC ${globalarrays_SOURCE_DIR}/ma)
target_include_directories(ga_sim PUBLIC ${globalarrays_SOURCE_DIR}/ga++/src)
target_include_directories(ga_sim PUBLIC ${globalarrays_SOURCE_DIR}/global/src)
target_include_directories(ga_sim PUBLIC ${globalarrays_BINARY_DIR}/ma)
target_include_directories(ga_sim PUBLIC ${globalarrays_BINARY_DIR}/gaf2c)
target_include_directories(ga_sim PUBLIC ${mxx_SOURCE_DIR}/include)

target_include_directories(ga_sim_no_tree PUBLIC ${globalarrays_SOURCE_DIR}/ma)
target_include_directories(ga_sim_no_tree PUBLIC ${globalarrays_SOURCE_DIR}/ga++/src)
target_include_directories(ga_sim_no_tree PUBLIC ${globalarrays_SOURCE_DIR}/global/src)
target_include_directories(ga_sim_no_tree PUBLIC ${globalarrays_BINARY_DIR}/ma)
target_include_directories(ga_sim_no_tree PUBLIC ${globalarrays_BINARY_DIR}/gaf2c)
target_include_directories(ga_sim_no_tree PUBLIC ${mxx_SOURCE_DIR}/include)

target_include_directories(ga_sim1D PUBLIC ${globalarrays_SOURCE_DIR}/ma)
target_include_directories(ga_sim1D PUBLIC ${globalarrays_SOURCE_DIR}/ga++/src)
target_include_directories(ga_sim1D PUBLIC ${globalarrays_SOURCE_DIR}/global/src)
target_include_directories(ga_sim1D PUBLIC ${globalarrays_BINARY_DIR}/ma)
target_include_directories(ga_sim1D PUBLIC ${globalarrays_BINARY_DIR}/gaf2c)

target_include_directories(ga_sim1D_hydro PUBLIC ${globalarrays_SOURCE_DIR}/ma)
target_include_directories(ga_sim1D_hydro PUBLIC ${globalarrays_SOURCE_DIR}/ga++/src)
target_include_directories(ga_sim1D_hydro PUBLIC ${globalarrays_SOURCE_DIR}/global/src)
target_include_directories(ga_sim1D_hydro PUBLIC ${globalarrays_BINARY_DIR}/ma)
target_include_directories(ga_sim1D_hydro PUBLIC ${globalarrays_BINARY_DIR}/gaf2c)

target_link_libraries(ga_sim GlobalArrays::ga MPI::MPI_CXX)
target_link_libraries(ga_sim_no_tree GlobalArrays::ga MPI::MPI_CXX)
target_link_libraries(ga_sim1D GlobalArrays::ga MPI::MPI_CXX)
target_link_libraries(ga_sim1D_hydro GlobalArrays::ga MPI::MPI_CXX)
target_link_libraries(omp_sim PUBLIC OpenMP::OpenMP_CXX)

