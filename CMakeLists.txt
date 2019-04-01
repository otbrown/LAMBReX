cmake_minimum_required(VERSION 3.5)

project(
  LAMBReX VERSION 0.1.0
  DESCRIPTION "Lattice Boltzmann code built on AMReX"
  LANGUAGES CXX
  )

option(LAMBREX_USE_CONAN "Build dependencies with Conan" ON)
option(LAMBREX_BUILD_EXAMPLES "Build the tests" ON)
option(LAMBREX_BUILD_TESTS "Build the tests" ON)

set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

if (LAMBREX_USE_CONAN)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup(TARGETS)
endif()

# You may need to set AMREX_DIR to help cmake find the config files
find_package(AMReX 19.03 REQUIRED)
if (NOT "${AMREX_DIM}" EQUAL 3)
  message(FATAL_ERROR "AMReX compiled not for 3D!")
endif()

add_subdirectory(src)
if (LAMBREX_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if (LAMBREX_BUILD_TESTS)
  # testing setup
  # Catch2 Test Framework
  find_package(Catch2 REQUIRED)
  enable_testing()
  include(CTest)
  include(Catch)

  add_subdirectory(tests)
endif()
