cmake_minimum_required(VERSION 3.11) # FetchContent is available in 3.11+
project(mapGen)
set(CMAKE_CXX_STANDARD 20)

# SFML 2.5
#find_package(SFML 2.5 COMPONENTS graphics audio window system network QUIET)

# Build types

if (UNIX)
  add_compile_options(-Wall -Wextra -Og)
  SET(CMAKE_BUILD_TYPE DEBUG)
  #Install sfml if not already installed
  if (NOT SFML_FOUND)
    include(FetchContent)

    FetchContent_Declare(
      SFML
      URL https://www.sfml-dev.org/files/SFML-2.5.1-linux-gcc-64-bit.tar.gz
    )
    FetchContent_GetProperties(SFML)
    if (NOT sfml_POPULATED)
      set(FETCHCONTENT_QUIET NO)
      FetchContent_Populate(sfml)
      set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
      include_directories("${sfml_SOURCE_DIR}/include")
    endif()
  endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()


# File to compile

add_executable(${PROJECT_NAME}
  osmiumGen.cpp
)


if (UNIX)
  target_link_libraries(${PROJECT_NAME}
    ${sfml_SOURCE_DIR}/lib/libsfml-graphics.so
    ${sfml_SOURCE_DIR}/lib/libsfml-audio.so
    ${sfml_SOURCE_DIR}/lib/libsfml-window.so
    ${sfml_SOURCE_DIR}/lib/libsfml-system.so
    ${sfml_SOURCE_DIR}/lib/libsfml-network.so
  )
endif()

# target_link_libraries(${PROJECT_NAME} sfml-graphics sfml-audio sfml-window sfml-system sfml-network)
set(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)


