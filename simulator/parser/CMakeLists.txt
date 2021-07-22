add_library(sim_parser STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/parser.cxx)
target_include_directories(sim_parser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(sim_parser PUBLIC biosoup::biosoup bioparser::bioparser)
set_target_default_warnings(sim_parser PRIVATE FALSE)
