cmake_minimum_required(VERSION 3.5)

project(
  LAMBReX VERSION 0.1.0
  DESCRIPTION "Lattice Boltzmann code built on AMReX"
  LANGUAGES CXX
  )
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LAMBREX_USE_CONAN "Build dependencies with Conan" OFF)
option(LAMBREX_BUILD_EXAMPLES "Build the examples" ON)
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
find_package(AMReX 19.08 REQUIRED 3D)

add_subdirectory(src)
if (LAMBREX_BUILD_EXAMPLES)
  add_subdirectory(examples)
  add_subdirectory(examples/debug)
endif()

if (LAMBREX_BUILD_TESTS)
  # testing setup
  # Catch2 Test Framework
  find_package(Catch2 REQUIRED)
  enable_testing()
  include(CTest)
  include(Catch)
  include(ParseAndAddCatchTests)

  add_subdirectory(tests)
endif()
