#include <catch2/catch.hpp>
#include <array>
#include "lambrex.h"

TEST_CASE("Single Value Initialisation", "[initialisation]")
{
  const double TEST_DENSITY = 0.63;
  const double TEST_VELOCITY = 0.23;
  const int TEST_NX = 11;
  const int TEST_NY = 12;
  const int TEST_NZ = 13;
  const int LEVEL = 0;
  std::array<int,NDIMS> periodicity({1,1,1});
  amrex::RealBox domain(AMREX_D_DECL(0.0,0.0,0.0), AMREX_D_DECL(1.0,1.0,1.0));

  AmrSim sim(TEST_NX, TEST_NY, TEST_NZ, 0.01, 0.01, periodicity, domain);
  sim.SetInitialDensity(TEST_DENSITY);
  sim.SetInitialVelocity(TEST_VELOCITY);
  sim.InitFromScratch(0.0);

  SECTION("check dimensions") {
    std::array<int,NDIMS> dims = sim.GetDims();
    REQUIRE(dims[0] == TEST_NX);
    REQUIRE(dims[1] == TEST_NY);
    REQUIRE(dims[2] == TEST_NZ);
  }

  SECTION("check density values") {
    for (int k = 0; k < TEST_NZ; ++k) {
      for (int j = 0; j < TEST_NY; ++j) {
        for (int i = 0; i < TEST_NX; ++i) {
          REQUIRE(sim.GetDensity(i, j, k, LEVEL) == Approx(TEST_DENSITY));
        }
      }
    }
  }

  SECTION("check velocity values") {
    for (int k = 0; k < TEST_NZ; ++k) {
      for (int j = 0; j < TEST_NY; ++j) {
        for (int i = 0; i < TEST_NX; ++i) {
          for (int dim = 0; dim < NDIMS; ++dim) {
            REQUIRE(sim.GetVelocity(i, j, k, dim, LEVEL) == Approx(TEST_VELOCITY));
          }
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
  const int NUMEL = TEST_NX * TEST_NY * TEST_NZ;
  const int LEVEL = 0;
  std::array<int,NDIMS> periodicity({1,1,1});
  std::vector<double> init_density(NUMEL);
  std::vector<double> init_velocity(3*NUMEL);
  amrex::RealBox domain(AMREX_D_DECL(0.0,0.0,0.0), AMREX_D_DECL(1.0,1.0,1.0));
  double val;

  val = 0.0;
  for (double& element : init_density) {
    element = val++;
  }

  val = 0.0;
  for (double& element : init_velocity) {
    element = val++;
  }

  AmrSim sim(TEST_NX, TEST_NY, TEST_NZ, 0.01, 0.01, periodicity, domain);
  sim.SetInitialDensity(init_density);
  sim.SetInitialVelocity(init_velocity);
  sim.InitFromScratch(0.0);

  SECTION("check dimensions") {
    std::array<int,NDIMS> dims = sim.GetDims();
    REQUIRE(dims[0] == TEST_NX);
    REQUIRE(dims[1] == TEST_NY);
    REQUIRE(dims[2] == TEST_NZ);
  }

  SECTION("check density values") {
    val = 0.0;
    // need C-like indexing here
    for (int i = 0; i < TEST_NX; ++i) {
      for (int j = 0; j < TEST_NY; ++j) {
        for (int k = 0; k < TEST_NZ; ++k) {
          REQUIRE(sim.GetDensity(i, j, k, LEVEL) == Approx(val++));
        }
      }
    }
  }

  SECTION("check velocity values") {
    val = 0.0;
    // need C-like indexing here
    for (int i = 0; i < TEST_NX; ++i) {
      for (int j = 0; j < TEST_NY; ++j) {
        for (int k = 0; k < TEST_NZ; ++k) {
          for (int dim = 0; dim < NDIMS; ++dim) {
            REQUIRE(sim.GetVelocity(i, j, k, dim, LEVEL) == Approx(val++));
          }
        }
      }
    }
  }
}
