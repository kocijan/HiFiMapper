include(FetchContent)

macro(get_remote_package NAME UPSTREAM TAG)
  FetchContent_Declare(
    ${NAME}
    GIT_REPOSITORY ${UPSTREAM}
    GIT_TAG ${TAG}
  )

  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.14)
    FetchContent_MakeAvailable(${NAME})
  else()
    FetchContent_GetProperties(${NAME}) 
    if (NOT ${NAME}_POPULATED)
      FetchContent_Populate(${NAME})
      add_subdirectory(
        ${${NAME}_SOURCE_DIR} 
        ${${NAME}_BINARY_DIR}
        EXCLUDE_FROM_ALL
      )
    endif()
  endif()
endmacro()
