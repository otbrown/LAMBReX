#include <catch2/catch.hpp>
#include <limits>
#include "lambrex.h"
#include "AmrTest.h"
#include <cstdio>

TEST_CASE("OneLevel", "[AMR]")
{
  const int NX = 10;
  const int NY = 10;
  const int NZ = 50;
  const int LEVEL = 0;
  const double DENSITY = 0.5;
  const double VELOCITY = 0.2;
  const double TAU = 0.01;

  lambrexSetAmr(NX, NY, NZ, LEVEL);
  AmrTest sim(TAU, TAU);
  sim.SetInitialDensity(DENSITY);
  sim.SetInitialVelocity(VELOCITY);
  sim.InitFromScratch(0.0);

  SECTION("ErrorEst") {
    const int LEVEL = 0;
    amrex::IntVect pos, lo, hi;
    amrex::TagBoxArray tba(sim.boxArray(LEVEL), sim.DistributionMap(LEVEL));
    tba.setVal(sim.boxArray(LEVEL), amrex::TagBox::SET);
    sim.CallErrorEst(LEVEL, tba);

    // at the moment ErrorEst should always tag every cell as clear (as we don't
    // have any control logic for when a cell should be refined..) so check that
    // this has happened
    for (amrex::MFIter mfi(tba); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      amrex::TagBox& tagfab = tba[mfi];

      lo = box.smallEnd();
      hi = box.bigEnd();
      for (int k = lo[2]; k <= hi[2]; ++k) {
        pos.setVal(2, k);
        for (int j = lo[1]; j <= hi[1]; ++j) {
          pos.setVal(1, j);
          for (int i = lo[0]; i <= hi[0]; ++i) {
            pos.setVal(0, i);
            INFO("(i j k) = (" << i << " " << j << " " << k << ")");
            REQUIRE(tagfab(pos) == amrex::TagBox::CLEAR);
          }
        }
      }
    }
  }

  SECTION("ClearLevel") {
    sim.CallClearLevel(LEVEL);

    REQUIRE(sim.DensityEmpty(LEVEL));
    REQUIRE(sim.VelocityEmpty(LEVEL));
    REQUIRE(sim.DistFnEmpty(LEVEL));
    REQUIRE(sim.GetTime(LEVEL) == Approx(0.0));
    REQUIRE(sim.GetTimeStep(LEVEL) == 0);
  }

  SECTION("MakeNewLevelFromScratch") {
    const double TIME = 1.0;
    const amrex::BoxArray ba(sim.boxArray(LEVEL));
    const amrex::DistributionMapping dm(sim.DistributionMap(LEVEL));

    sim.CallClearLevel(LEVEL);
    sim.CallMakeNewLevelFromScratch(ba, dm, TIME);

    REQUIRE_FALSE(sim.DensityEmpty(LEVEL));
    REQUIRE_FALSE(sim.VelocityEmpty(LEVEL));
    REQUIRE_FALSE(sim.DistFnEmpty(LEVEL));

    REQUIRE(sim.GetTime(LEVEL) == TIME);

    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
          REQUIRE(sim.GetDensity(i, j, k, LEVEL) == Approx(DENSITY));
          for (int dim = 0; dim < NDIMS; ++dim)
            REQUIRE(sim.GetVelocity(i, j, k, dim, LEVEL) == Approx(VELOCITY));
        }
      }
    }
  }

  SECTION("RemakeLevel") {
    const double TIME = sim.GetTime(LEVEL);
    const amrex::BoxArray ba(sim.boxArray(LEVEL));
    const amrex::DistributionMapping dm(sim.DistributionMap(LEVEL));

    sim.CallRemakeLevel(LEVEL, TIME, ba, dm);

    // This is a tougher test than it appears, as density and velocity are
    // recalculated from the distribution function
    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
          REQUIRE(sim.GetDensity(i, j, k, LEVEL) == Approx(DENSITY));
          for (int dim = 0; dim < NDIMS; ++dim)
            REQUIRE(sim.GetVelocity(i, j, k, dim, LEVEL) == Approx(VELOCITY));
        }
      }
    }
  }
}

TEST_CASE("TwoLevel", "[AMR]") {
  const int NX = 48;
  const int NY = 24;
  const int NZ = 12;
  const int MAX_LEVEL = 1;
  const int NUM_LEVELS = MAX_LEVEL + 1;
  const double DENSITY = 0.8;
  const double VELOCITY = 0.1;
  const double TAU = 0.01;

  lambrexSetAmr(NX, NY, NZ, MAX_LEVEL);
  AmrTest sim(TAU, TAU);

  sim.SetInitialDensity(DENSITY);
  sim.SetInitialVelocity(VELOCITY);
  sim.InitFromScratch(0.0);

  SECTION("Initialisation") {
    REQUIRE(sim.maxLevel() == MAX_LEVEL);
    REQUIRE(sim.GetDensity().size() == NUM_LEVELS);
    REQUIRE(sim.GetVelocity().size() == NUM_LEVELS);
    REQUIRE(sim.GetDistFn().size() == NUM_LEVELS);
    REQUIRE(sim.GetSimTime().size() == NUM_LEVELS);
    REQUIRE(sim.GetDt().size() == NUM_LEVELS);
    REQUIRE(sim.GetTimeStep().size() == NUM_LEVELS);

    // refinement ratio size is NUM_LEVELS-1, since only has meaning "between"
    // levels. Default is 2 in every direction, we don't change this yet.
    REQUIRE(sim.refRatio().size() == MAX_LEVEL);
    for (int level = 1; level < NUM_LEVELS; ++level) {
      for (int dim = 0; dim < NDIMS; ++dim)
        REQUIRE(sim.refRatio(level-1)[dim] == 2);
    }

    REQUIRE_FALSE(sim.DensityEmpty(0));
    REQUIRE_FALSE(sim.VelocityEmpty(0));
    REQUIRE_FALSE(sim.DistFnEmpty(0));

    REQUIRE(sim.DensityEmpty(MAX_LEVEL));
    REQUIRE(sim.VelocityEmpty(MAX_LEVEL));
    REQUIRE(sim.DistFnEmpty(MAX_LEVEL));
  }
}
