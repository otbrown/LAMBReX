#include <catch2/catch.hpp>
#include <limits>
#include "lambrex.h"
#include "pulseRegression.h"

TEST_CASE("pulse Regression", "[regression]")
{
  const int NX = 10;
  const int NY = 10;
  const int NZ = 50;
  const double TAU = 0.5;
  double amplitude = 0.01;
  int periodicity[3] = {1, 1, 1};
  double z_mean, tmp;
  int i, j, k, n, rhodex, veldex;

  Simulation sim(NX, NY, NZ, TAU, TAU, periodicity);

  // set initial density
  sim.setDensity(1.0);

  // set enhanced density
  k = NZ/2 - 1;
  for (j = 0; j < NY; ++j) {
    for (i = 0; i < NX; ++i) {
      sim.setDensity(i, j, k, 1.0 + amplitude);
    }
  }

  // renormalise
  z_mean = 0.0;
  for (k = 0; k < NZ; ++k) z_mean += sim.getDensity(NX/2, NY/2, k);
  z_mean /= NZ;
  for (k = 0; k < NZ; ++k) {
    for (j = 0; j < NY; ++j) {
      for (i = 0; i < NX; ++i) {
        tmp = sim.getDensity(i, j, k);
        sim.setDensity(i, j, k, tmp / z_mean);
      }
    }
  }

  // set initial velocity
  sim.setVelocity(0.0);

  // calculate initial distribution function
  sim.calcEquilibriumDist();

  rhodex = 0;
  veldex = 0;
  for (k = 0; k < NZ; ++k) {
    for (j = 0; j < NY; ++j) {
      for (i = 0; i < NX; ++i) {
        INFO("t=0 (i,j,k)=(" << i << "," << j << "," << k << ") rhodex="
             << rhodex)
        REQUIRE(sim.getDensity(i, j, k) == Approx(RHO_t0[rhodex++]));
        for (n = 0; n < AMREX_SPACEDIM; ++n) {
          INFO("veldex=" << veldex << " n=" << n)
          REQUIRE(sim.getVelocity(i, j, k, n) == Approx(VEL_t0[veldex++]));
        }
      }
    }
  }

  sim.iterate(100);
  sim.calcHydroVars();
  rhodex = 0;
  veldex = 0;
  for (k = 0; k < NZ; ++k) {
    for (j = 0; j < NY; ++j) {
      for (i = 0; i < NX; ++i) {
        INFO("t=100 (i,j,k)=(" << i << "," << j << "," << k << ") rhodex="
             << rhodex)
        REQUIRE(sim.getDensity(i, j, k) == Approx(RHO_t100[rhodex++]));
        for (n = 0; n < AMREX_SPACEDIM; ++n) {
          INFO("veldex=" << veldex << " n=" << n)
          REQUIRE(sim.getVelocity(i, j, k, n) == Approx(VEL_t100[veldex++]));
        }
      }
    }
  }

  sim.iterate(100);
  sim.calcHydroVars();
  rhodex = 0;
  veldex = 0;
  for (k = 0; k < NZ; ++k) {
    for (j = 0; j < NY; ++j) {
      for (i = 0; i < NX; ++i) {
        INFO("t=200 (i,j,k)=(" << i << "," << j << "," << k << ") rhodex="
             << rhodex)
        REQUIRE(sim.getDensity(i, j, k) == Approx(RHO_t200[rhodex++]));
        for (n = 0; n < AMREX_SPACEDIM; ++n) {
          INFO("veldex=" << veldex << " n=" << n)
          REQUIRE(sim.getVelocity(i, j, k, n) == Approx(VEL_t200[veldex++]));
        }
      }
    }
  }
}
