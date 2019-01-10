 # Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.0.2)

enable_language(Fortran)

# Project's name
project(davidson_eigensolver)

set( CMAKE_VERBOSE_MAKEFILE on )

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})


# search for LAPACK
find_package(LAPACK COMPONENTS Fortran REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})
MESSAGE (STATUS "LAPACK PATH:" ${LAPACK_LIBRARIES})
set(LINEAR_ALGEBRA ${LAPACK_LIBRARIES})

# # search for Doxygen
# find_package(DOXYGEN COMPONENTS REQUIRED)
# MESSAGE(STATUS "DOXYGEN VERSION: " ${DOXYGEN_VERSION})

add_subdirectory(src)
