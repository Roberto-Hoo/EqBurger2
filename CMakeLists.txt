cmake_minimum_required(VERSION 3.19)
project(EqBurger2)

set(CMAKE_CXX_STANDARD 14)

add_executable(EqBurger2 EqBurger2.cpp)
# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpic++)
set(SOURCE_FILES EqBurger2.cpp)
