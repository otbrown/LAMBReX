cmake_minimum_required(VERSION 3.13)

project(
  LAMBReX VERSION 0.1.0
  DESCRIPTION "Lattice Boltzmann code built on AMReX"
  LANGUAGES CXX
)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

# AMReX do not provide FindAMReX.cmake (yet)
set(AMReX_DIR "${CMAKE_SOURCE_DIR}/libamrex/lib/cmake/AMReX")
find_package(AMReX 19.02 REQUIRED)

# we could get this from user, but for now hardcode as MUST be 3
set(DIM "3")

add_subdirectory(src)
add_subdirectory(examples)

# testing setup
enable_testing()

add_test(
  NAME "pulse_diff"
  COMMAND ${CMAKE_COMMAND} -E compare_files pulse.out blessed_pulse.out
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bins
)
