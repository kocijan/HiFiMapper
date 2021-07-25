cmake_minimum_required(VERSION 3.11)
project(HifFiMapper VERSION 1.0.1
                    LANGUAGES C CXX
                    DESCRIPTION "HifFiMapper is a reference assembly program optimized for PacBio-Hifi reads")

include(FetchContent)
include(GNUInstallDirs)

include(${PROJECT_SOURCE_DIR}/cmake/Content.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/Warnings.cmake)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(bioparser 3.0.13 QUIET)
if (NOT bioparser_FOUND)
  get_remote_package(
    bioparser 
    https://github.com/rvaser/bioparser
    3.0.13
  )
endif()

find_package(biosoup 0.10.0 QUIET)
if (NOT biosoup_FOUND)
  get_remote_package(
    biosoup  
    https://github.com/rvaser/biosoup
    0.10.0
  )
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -O0 -g")
  set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -O0 -g")
  message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mapper)
if (BUILD_SIMULATOR)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/simulator)
endif()

