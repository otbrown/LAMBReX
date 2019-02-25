#include <catch2/catch.hpp>
#include <array>
#include <cstdio>
#include "lambrex.h"

TEST_CASE("Single Value Initialisation", "[initialisation]")
{
  const double TEST_DENSITY = 0.63;
  const int TEST_NX = 11;
  const int TEST_NY = 12;
  const int TEST_NZ = 13;
  int periodicity[3] = {1,1,1};

  Simulation sim(TEST_NX, TEST_NY, TEST_NZ, 0.01, 0.01, periodicity);
  sim.setDensity(TEST_DENSITY);

  SECTION("check dimensions") {
    std::array<int,3> dims = sim.getDims();
    REQUIRE(dims[0] == TEST_NX);
    REQUIRE(dims[1] == TEST_NY);
    REQUIRE(dims[2] == TEST_NZ);
  }

  SECTION("check density values") {
    for (int k = 0; k < TEST_NZ; ++k) {
      for (int j = 0; j < TEST_NY; ++j) {
        for (int i = 0; i < TEST_NX; ++i) {
          REQUIRE(sim.getDensity(i, j, k) == Approx(TEST_DENSITY));
        }
      }
    }
  }
}

TEST_CASE("Element-wise Initialisation", "[initialisation]")
{
  const int TEST_NX = 16;
  const int TEST_NY = 9;
  const int TEST_NZ = 8;
  int periodicity[3] = {1,1,1};

  Simulation sim(TEST_NX, TEST_NY, TEST_NZ, 0.01, 0.01, periodicity);

  double val = 0.0;
  for (int k = 0; k < TEST_NZ; ++k) {
    for (int j = 0; j < TEST_NY; ++j) {
      for (int i = 0; i < TEST_NX; ++i) {
        sim.setDensity(i, j, k, val++);
      }
    }
  }

  SECTION("check dimensions") {
    std::array<int,3> dims = sim.getDims();
    REQUIRE(dims[0] == TEST_NX);
    REQUIRE(dims[1] == TEST_NY);
    REQUIRE(dims[2] == TEST_NZ);
  }

  SECTION("check density values") {
    val = 0.0;
    for (int k = 0; k < TEST_NZ; ++k) {
      for (int j = 0; j < TEST_NY; ++j) {
        for (int i = 0; i < TEST_NX; ++i) {
          REQUIRE(sim.getDensity(i, j, k) == Approx(val++));
        }
      }
    }
  }
}
