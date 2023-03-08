cmake_minimum_required(VERSION 3.24)
project(bcg_method)

find_package(CBLAS)

set(CMAKE_CXX_STANDARD 17)
set(GCC_COVERAGE_COMPILE_FLAGS "-lcblas")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(BLAS REQUIRED)


add_executable(bcg_method main.cpp)
