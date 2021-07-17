# CMake configuration 

cmake_minimum_required(VERSION 3.8)

project(HiFi-Simulator VERSION 1.0)

# Required C++ standard.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add the executable
add_executable(Refgen reference_generator.cpp)
add_executable(Readsgen reads_generator.cpp)

add_subdirectory (../vendor/bioparser)
add_subdirectory (../vendor/biosoup)

if (NOT TARGET bioparser)
  add_subdirectory(../vendor/bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(Readsgen bioparser::bioparser)

if (NOT TARGET biosoup)
  add_subdirectory(../vendor/biosoup EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(Readsgen biosoup::biosoup)
