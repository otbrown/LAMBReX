#include <catch2/catch.hpp>
#include <limits>
#include "lambrex.h"
#include "pulseRegression.h"

TEST_CASE("pulse Regression", "[regression]")
{
  const int NX = 10;
  const int NY = 10;
  const int NZ = 50;
  const int NUMEL = NX*NY*NZ;
  const double TAU = 0.5;
  const double AMPLITUDE = 0.01;
  double z_mean;
  int i, j, k, n, rhodex, veldex;
  const int LEVEL = 0;

  // initial density
  std::vector<double> rho_0(NUMEL, 1.0);
  k = NZ / 2;
  for (i = 0; i < NX; ++i) {
    for (j = 0; j < NY; ++j) {
      rho_0[i*NY*NZ + j*NZ + k - 1] += AMPLITUDE;
    }
  }

  z_mean = 0.0;
  for (k = 0; k < NZ; ++k) z_mean += rho_0[(NUMEL + NY*NZ) / 2 + k];
  z_mean /= NZ;
  for (i = 0; i < NUMEL; ++i) rho_0[i] /= z_mean;

  // initial velocity
  double u_0 = 0.0;

  // initialise AmrSim
  lambrexSetAmr(NX, NY, NZ, LEVEL);
  AmrSim sim(TAU, TAU);
  sim.SetInitialDensity(rho_0);
  sim.SetInitialVelocity(u_0);
  sim.InitFromScratch(0.0);

  rhodex = 0;
  veldex = 0;
  for (k = 0; k < NZ; ++k) {
    for (j = 0; j < NY; ++j) {
      for (i = 0; i < NX; ++i) {
        INFO("t=0 (i,j,k)=(" << i << "," << j << "," << k << ") rhodex="
             << rhodex)
        REQUIRE(sim.GetDensity(i, j, k, LEVEL) == Approx(RHO_t0[rhodex++]));
        for (n = 0; n < AMREX_SPACEDIM; ++n) {
          INFO("veldex=" << veldex << " n=" << n)
          REQUIRE(sim.GetVelocity(i, j, k, n, LEVEL) == Approx(VEL_t0[veldex++]));
        }
      }
    }
  }

  sim.Iterate(100);
  sim.CalcHydroVars(LEVEL);
  rhodex = 0;
  veldex = 0;
  for (k = 0; k < NZ; ++k) {
    for (j = 0; j < NY; ++j) {
      for (i = 0; i < NX; ++i) {
        INFO("t=100 (i,j,k)=(" << i << "," << j << "," << k << ") rhodex="
             << rhodex)
        REQUIRE(sim.GetDensity(i, j, k, LEVEL) == Approx(RHO_t100[rhodex++]));
        for (n = 0; n < AMREX_SPACEDIM; ++n) {
          INFO("veldex=" << veldex << " n=" << n)
          REQUIRE(sim.GetVelocity(i, j, k, n, LEVEL) == Approx(VEL_t100[veldex++]));
        }
      }
    }
  }

  sim.Iterate(100);
  sim.CalcHydroVars(LEVEL);
  rhodex = 0;
  veldex = 0;
  for (k = 0; k < NZ; ++k) {
    for (j = 0; j < NY; ++j) {
      for (i = 0; i < NX; ++i) {
        INFO("t=200 (i,j,k)=(" << i << "," << j << "," << k << ") rhodex="
             << rhodex)
        REQUIRE(sim.GetDensity(i, j, k, LEVEL) == Approx(RHO_t200[rhodex++]));
        for (n = 0; n < AMREX_SPACEDIM; ++n) {
          INFO("veldex=" << veldex << " n=" << n)
          REQUIRE(sim.GetVelocity(i, j, k, n, LEVEL) == Approx(VEL_t200[veldex++]));
        }
      }
    }
  }
}

TEST_CASE("ml_pulse Regression", "[regression]")
{
  const int NX = 10;
  const int NY = 10;
  const int NZ = 50;
  const int NUMEL = NX*NY*NZ;
  const double TAU = 0.5;
  const int MAX_LEVEL = 1;
  const double AMPLITUDE = 0.01;
  const std::array<int,3> LO = { 0, 0, 0 };
  const std::array<int,3> HI = { NX-1, NY-1, NZ-1 };
  double z_mean;
  int i, j, k, n, rhodex, veldex, level, ratio;

  // initial density
  std::vector<double> rho_0(NUMEL, 1.0);
  k = NZ / 2;
  for (i = 0; i < NX; ++i) {
    for (j = 0; j < NY; ++j) {
      rho_0[i*NY*NZ + j*NZ + k - 1] += AMPLITUDE;
    }
  }

  z_mean = 0.0;
  for (k = 0; k < NZ; ++k) z_mean += rho_0[(NUMEL + NY*NZ) / 2 + k];
  z_mean /= NZ;
  for (i = 0; i < NUMEL; ++i) rho_0[i] /= z_mean;

  // initial velocity
  double u_0 = 0.0;

  // initialise AmrSim
  lambrexSetAmr(NX, NY, NZ, MAX_LEVEL);
  AmrSim sim(TAU, TAU);
  sim.SetInitialDensity(rho_0);
  sim.SetInitialVelocity(u_0);
  sim.InitFromScratch(0.0);
  sim.SetStaticRefinement(0, LO, HI);

  ratio = 1;
  for (level = 0; level <= MAX_LEVEL; ++level) {
    rhodex = 0;
    veldex = 0;
    for (k = ratio*LO[2]; k <= ratio*HI[2]; k += ratio) {
      for (j = ratio*LO[1]; j <= ratio*HI[1]; j += ratio) {
        for (i = ratio*LO[0]; i <= ratio*HI[0]; i += ratio) {
          INFO("t=0 (i,j,k)=(" << i << "," << j << "," << k << ") rhodex="
               << rhodex)
          REQUIRE(sim.GetDensity(i, j, k, level) == Approx(RHO_t0[rhodex++]));
          for (n = 0; n < AMREX_SPACEDIM; ++n) {
            INFO("veldex=" << veldex << " n=" << n)
            REQUIRE(sim.GetVelocity(i, j, k, n, level) == Approx(VEL_t0[veldex++]));
          }
        }
      }
    }
    ratio *= 2;
  }

  level = 0;
  sim.Iterate(100);
  sim.CalcHydroVars(level);

  rhodex = 0;
  veldex = 0;
  for (k = LO[2]; k <= HI[2]; ++k) {
    for (j = LO[1]; j <= HI[1]; ++j) {
      for (i = LO[0]; i <= HI[0]; ++i) {
        INFO("t=100 (i,j,k)=(" << i << "," << j << "," << k << ") level="
          << level << " rhodex=" << rhodex)
        REQUIRE(sim.GetDensity(i, j, k, level) == Approx(RHO_t100[rhodex++]));
        for (n = 0; n < AMREX_SPACEDIM; ++n) {
          INFO("veldex=" << veldex << " n=" << n)
          REQUIRE(sim.GetVelocity(i, j, k, n, level) == Approx(VEL_t100[veldex++]));
        }
      }
    }
  }

  sim.Iterate(100);
  for (level = 0; level <= MAX_LEVEL; ++level) sim.CalcHydroVars(level);

  ratio = 1;
  for (level = 0; level <= MAX_LEVEL; ++level) {
    rhodex = 0;
    veldex = 0;
    for (k = ratio*LO[2]; k <= ratio*HI[2]; k += ratio) {
      for (j = ratio*LO[1]; j <= ratio*HI[1]; j += ratio) {
        for (i = ratio*LO[0]; i <= ratio*HI[0]; i += ratio) {
          INFO("t=100 (i,j,k)=(" << i << "," << j << "," << k << ") rhodex="
               << rhodex)
          REQUIRE(sim.GetDensity(i, j, k, level) == Approx(RHO_t200[rhodex++]));
          for (n = 0; n < AMREX_SPACEDIM; ++n) {
            INFO("veldex=" << veldex << " n=" << n)
            REQUIRE(sim.GetVelocity(i, j, k, n, level) == Approx(VEL_t200[veldex++]));
          }
        }
      }
    }
    ratio *= 2;
  }
}
