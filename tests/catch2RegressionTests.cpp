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
  double amplitude = 0.01;
  std::array<int,3> periodicity = {1, 1, 1};
  double z_mean;
  int i, j, k, n, rhodex, veldex;
  const int LEVEL = 0;
  amrex::RealBox domain(AMREX_D_DECL(0.0,0.0,0.0), AMREX_D_DECL(1.0,1.0,1.0));

  // initial density
  std::vector<double> rho_0(NUMEL, 1.0);
  k = NZ / 2;
  for (i = 0; i < NX; ++i) {
    for (j = 0; j < NY; ++j) {
      rho_0[i*NY*NZ + j*NZ + k - 1] += amplitude;
    }
  }

  z_mean = 0.0;
  for (k = 0; k < NZ; ++k) z_mean += rho_0[(NUMEL + NY*NZ) / 2 + k];
  z_mean /= NZ;
  for (i = 0; i < NUMEL; ++i) rho_0[i] /= z_mean;

  // initial velocity
  double u_0 = 0.0;

  // initialise AmrSim
  AmrSim sim(NX, NY, NZ, TAU, TAU, periodicity, domain);
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
