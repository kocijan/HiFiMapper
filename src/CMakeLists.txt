# CMake configuration 

cmake_minimum_required(VERSION 3.8)

project(HiFi-Mapper VERSION 1.0)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(RELEASE_FLAGS "-Ofast -DNDEBUG -DBOOST_DISABLE_ASSERTS")
set(DEBUG_FLAGS "-O0 -ggdb3")

set(CMAKE_CXX_FLAGS_RELEASE ${RELEASE_FLAGS})
set(CMAKE_C_FLAGS_RELEASE ${RELEASE_FLAGS})

set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_C_FLAGS_DEBUG ${DEBUG_FLAGS})


# Required C++ standard.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add the executable
add_executable(HiFiMapper main.cpp)

find_package(bioparser)
target_link_libraries(HiFiMapper bioparser::bioparser)

find_package(biosoup)
target_link_libraries(HiFiMapper biosoup::biosoup)
