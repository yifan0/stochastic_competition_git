cmake_minimum_required (VERSION 3.18.0 FATAL_ERROR)
project (stochastic_competition VERSION 1.0.0 LANGUAGES CXX)

find_package (MPI REQUIRED)

set(ENABLE_FORTRAN OFF)

include(FetchContent)
FetchContent_Declare(
    globalarrays
    GIT_REPOSITORY https://github.com/GlobalArrays/ga.git
    GIT_TAG "v5.8.2"
)
FetchContent_MakeAvailable(globalarrays)

include_directories(agner_random include)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++17")

set(targets
    sim
    tree_sim
    mpi_sim
    ga_sim
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

target_link_libraries(mpi_sim MPI::MPI_CXX)
target_link_libraries(ga_sim GlobalArrays::ga MPI::MPI_CXX)
