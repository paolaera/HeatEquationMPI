cmake_minimum_required(VERSION 3.14)
project(HeatEquation C)


find_library(MATLAB libmat.so PATHS /usr/local/MATLAB/R2019b/bin/glnxa64/)
find_library(MATLAB libmx.so PATHS /usr/local/MATLAB/R2019b/bin/glnxa64/)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

set(SOURCE_FILES main.c cmake-build-debug/MyLibrary.c cmake-build-debug/MyLibrary.h cmake-build-debug/MyLibrary.c cmake-build-debug/MyLibrary.h)
add_executable(HeatEquation ${SOURCE_FILES})


target_link_libraries(HeatEquation libmat.so "${MATLAB}")
target_link_libraries(HeatEquation libmx.so "${MATLAB}")
target_link_libraries(HeatEquation m)