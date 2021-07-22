find_package(OpenMP REQUIRED)

set(MAPPER_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mapper/chaining.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mapper/suffix_array.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mapper/parser.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mapper/lcp.cxx

  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cxx
)

add_executable(mapper ${MAPPER_SOURCES})
target_include_directories(mapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(mapper PRIVATE 
  biosoup::biosoup 
  bioparser::bioparser

  OpenMP::OpenMP_C
  OpenMP::OpenMP_CXX
)

set_target_default_warnings(mapper PRIVATE FALSE)


