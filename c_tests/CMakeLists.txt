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

add_subdirectory(agner_random)
include_directories(include)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++17")

set(targets
    sim
)

foreach(TARGET ${targets})
    add_executable(${TARGET} "${TARGET}.cpp")
    target_compile_features(${TARGET} PUBLIC cxx_std_17)
    target_include_directories(${TARGET} PUBLIC ${PROJECT_SOURCE_DIR}/agner_random)
    target_include_directories(${TARGET} PUBLIC ${globalarrays_SOURCE_DIR}/ma)
    target_include_directories(${TARGET} PUBLIC ${globalarrays_SOURCE_DIR}/ga++/src)
    target_include_directories(${TARGET} PUBLIC ${globalarrays_SOURCE_DIR}/global/src)
    target_include_directories(${TARGET} PUBLIC ${globalarrays_BINARY_DIR}/ma)
    target_include_directories(${TARGET} PUBLIC ${globalarrays_BINARY_DIR}/gaf2c)
    target_include_directories(${TARGET} PUBLIC ${mxx_SOURCE_DIR}/include)
    target_link_libraries(${TARGET} GlobalArrays::ga MPI::MPI_CXX Agner)
endforeach(TARGET)

