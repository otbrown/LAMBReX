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
    sim.CallMakeNewLevelFromScratch(LEVEL, ba, dm, TIME);

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
    const double TIME = 2.5;
    const amrex::BoxArray ba(sim.boxArray(LEVEL));
    const amrex::DistributionMapping dm(sim.DistributionMap(LEVEL));

    sim.CallRemakeLevel(LEVEL, TIME, ba, dm);

    // This is a tougher test than it appears, as density and velocity are
    // recalculated from the distribution function
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
  const std::array<int,NDIMS> LO_CORNER({0, 0, 0});
  const std::array<int,NDIMS> HI_CORNER({NX-1, NY-1, NZ-1});

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

  SECTION("MakeNewLevelFromCoarse") {
    // use box array and distribution mapping from level 0
    const amrex::BoxArray ba(sim.boxArray(0));
    const amrex::DistributionMapping dm(sim.DistributionMap(0));

    for (int level = 1; level <= MAX_LEVEL; ++level) {
      sim.CallMakeNewLevelFromCoarse(level, ba, dm);
    }

    for (int level = 0; level <= MAX_LEVEL; ++level)
      REQUIRE_FALSE(sim.DensityEmpty(level));

    for (int level = 0; level <= MAX_LEVEL; ++level) {
      REQUIRE(sim.GetTime(level) == 0.0);
      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {
            REQUIRE(sim.GetDensity(i, j, k, level) == Approx(DENSITY));
            for (int dim = 0; dim < NDIMS; ++dim)
              REQUIRE(sim.GetVelocity(i, j, k, dim, level) == Approx(VELOCITY));
          }
        }
      }
    }
  }

  SECTION("ClearLevel") {
    int level;

    // make sure there is something to clear on each level
    const amrex::BoxArray ba(sim.boxArray(0));
    const amrex::DistributionMapping dm(sim.DistributionMap(0));
    for (level = 1; level <= MAX_LEVEL; ++level)
      sim.CallMakeNewLevelFromCoarse(level, ba, dm);

    // now clear them
    for (level = 0; level <= MAX_LEVEL; ++level) sim.CallClearLevel(level);

    for (level = 0; level <= MAX_LEVEL; ++level) {
      REQUIRE(sim.DensityEmpty(level));
      REQUIRE(sim.VelocityEmpty(level));
      REQUIRE(sim.DistFnEmpty(level));
      REQUIRE(sim.GetTime(level) == 0.0);
      REQUIRE(sim.GetTimeStep(level) == 0);
    }
  }

  SECTION("ErrorEst") {
    int level;
    amrex::IntVect pos, lo, hi;
    const amrex::BoxArray ba(sim.boxArray(0));
    const amrex::DistributionMapping dm(sim.DistributionMap(0));
    std::vector<amrex::TagBoxArray> tba;

    for (level = 1; level <= MAX_LEVEL; ++level)
      sim.CallMakeNewLevelFromCoarse(level, ba, dm);

    for (level = 0; level <= MAX_LEVEL; ++level) {
      tba.emplace_back(ba, dm);
      tba.at(level).setVal(sim.boxArray(level), amrex::TagBox::SET);
      sim.CallErrorEst(level, tba.at(level));
    }

    for (level = 0; level <= MAX_LEVEL; ++level) {
      // at the moment ErrorEst should always tag every cell as clear
      for (amrex::MFIter mfi(tba.at(level)); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        amrex::TagBox& tagfab = tba.at(level)[mfi];

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

    // now set a static refinement area...
    for (level = 0; level < MAX_LEVEL; ++level)
      sim.SetStaticRefinement(level, LO_CORNER, HI_CORNER);

    // and call ErrorEst
    for (level = 0; level < MAX_LEVEL; ++level)
      sim.CallErrorEst(level, tba.at(level));

    // REQUIRE that all tags on all levels except the finest are now SET
    for (level = 0; level < MAX_LEVEL; ++level) {
      for (amrex::MFIter mfi(tba.at(level)); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        amrex::TagBox& tagfab = tba.at(level)[mfi];

        lo = box.smallEnd();
        hi = box.bigEnd();
        for (int k = lo[2]; k <= hi[2]; ++k) {
          pos.setVal(2, k);
          for (int j = lo[1]; j <= hi[1]; ++j) {
            pos.setVal(1, j);
            for (int i = lo[0]; i <= hi[0]; ++i) {
              pos.setVal(0, i);
              INFO("(i j k) = (" << i << " " << j << " " << k << ")");
              REQUIRE(tagfab(pos) == amrex::TagBox::SET);
            }
          }
        }
      }
    }

    // unset static refinements and ErrorEst again
    for (level = 0; level < MAX_LEVEL; ++level) {
      sim.UnsetStaticRefinement(level);
      sim.CallErrorEst(level, tba.at(level));
    }

    // REQUIRE that all tags on all levels except the finest are CLEAR again
    for (level = 0; level < MAX_LEVEL; ++level) {
      for (amrex::MFIter mfi(tba.at(level)); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        amrex::TagBox& tagfab = tba.at(level)[mfi];

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
  }

  SECTION("MakeNewLevelFromScratch") {
    int level;
    const double TIME = 1.0;
    const amrex::BoxArray ba(sim.boxArray(0));
    const amrex::DistributionMapping dm(sim.DistributionMap(0));

    // empty out sim first
    for (level = 0; level <= MAX_LEVEL; ++level) sim.CallClearLevel(level);

    // now make new levels
    for (level = 0; level <= MAX_LEVEL; ++level)
      sim.CallMakeNewLevelFromScratch(level, ba, dm, TIME);

    // coarsest level should have been fully initialised
    level = 0;
    REQUIRE_FALSE(sim.DensityEmpty(level));
    REQUIRE_FALSE(sim.VelocityEmpty(level));
    REQUIRE_FALSE(sim.DistFnEmpty(level));

    REQUIRE(sim.GetTime(level) == TIME);
    REQUIRE(sim.GetTimeStep(level) == 0);

    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
          REQUIRE(sim.GetDensity(i, j, k, level) == Approx(DENSITY));
          for (int dim = 0; dim < NDIMS; ++dim)
            REQUIRE(sim.GetVelocity(i, j, k, dim, level) == Approx(VELOCITY));
        }
      }
    }

    // all other levels should be not empty, but uninitialised
    for (level = 1; level <= MAX_LEVEL; ++level) {
      REQUIRE_FALSE(sim.DensityEmpty(level));
      REQUIRE_FALSE(sim.VelocityEmpty(level));
      REQUIRE_FALSE(sim.DistFnEmpty(level));
      REQUIRE(sim.GetTime(level) == TIME);
      REQUIRE(sim.GetTimeStep(level) == 0);
    }
  }

  SECTION("RemakeLevel") {
    int level;
    const double TIME = 1.3;
    const amrex::BoxArray ba(sim.boxArray(0));
    const amrex::DistributionMapping dm(sim.DistributionMap(0));

    // make sure all levels are initialised
    for (level = 1; level <= MAX_LEVEL; ++level)
      sim.CallMakeNewLevelFromCoarse(level, ba, dm);

    // remake all levels and REQUIRE
    for (level = 0; level <= MAX_LEVEL; ++level)
      sim.CallRemakeLevel(level, TIME, ba, dm);

    for (level = 0; level <= MAX_LEVEL; ++level) {
      REQUIRE(sim.GetTime(level) == TIME);
      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {
            REQUIRE(sim.GetDensity(i, j, k, level) == Approx(DENSITY));
            for (int dim = 0; dim < NDIMS; ++dim)
              REQUIRE(sim.GetVelocity(i, j, k, dim, level) == Approx(VELOCITY));
          }
        }
      }
    }
  }
}
