cmake_minimum_required(VERSION 3.13)

project(LAMBReX VERSION 0.1.0
                DESCRIPTION "Lattice Boltzmann code built on AMReX"
                LANGUAGES CXX)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

add_subdirectory(src)
add_subdirectory(examples)
