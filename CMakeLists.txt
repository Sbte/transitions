cmake_minimum_required(VERSION 3.5)
project(transitions CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS            "${CMAKE_CXX_FLAGS} -O3 -g")

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(main src/main.cpp src/timestepper.cpp)
target_link_libraries(main blas)
