add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/parser)

add_executable(Refgen ${CMAKE_CURRENT_SOURCE_DIR}/reference_generator.cpp)
add_executable(Readsgen ${CMAKE_CURRENT_SOURCE_DIR}/reads_generator.cpp)

target_link_libraries(Readsgen bioparser::bioparser biosoup::biosoup sim_parser)
