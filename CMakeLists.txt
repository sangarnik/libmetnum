cmake_minimum_required(VERSION 2.8.12)

project(metnum CXX C)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(lib/armadillo-12.2.0)

add_executable(metnum 
  src/main.cpp
  src/metnum_common.cpp
  src/metnum_diffeq.cpp
  src/metnum_eqsolver.cpp
  src/metnum_rootfinder.cpp
)

target_link_libraries(metnum armadillo)
target_include_directories(metnum PUBLIC include)
target_include_directories(metnum PUBLIC lib/armadillo-12.2.0/include)
