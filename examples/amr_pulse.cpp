#include <iostream>
#include "lambrex.h"

void printDensity(const AmrSim&);
void printVelocity(const AmrSim&);

int main (int argc, char * argv[])
{
  const int nx = 10;
  const int ny = 10;
  const int nz = 50;
  const double tau = 0.5;
  double amplitude = 0.01;
  std::array<int,3> periodicity = {1, 1, 1};
  const int numel = nx * ny * nz;
  int i, j, k;
  double z_mean;

  // initial density
  std::vector<double> rho_0(numel, 1.0);
  k = nz / 2;
  for (i = 0; i < nx; ++i) {
    for (j = 0; j < ny; ++j) {
      rho_0[i*ny*nz + j*nz + k - 1] += amplitude;
    }
  }

  z_mean = 0.0;
  for (k = 0; k < nz; ++k) z_mean += rho_0[(numel + ny*nz) / 2 + k];
  z_mean /= nz;
  for (i = 0; i < numel; ++i) rho_0[i] /= z_mean;

  // initial velocity
  double u_0 = 0.0;

  amrex::Initialize(argc, argv, false);
  {
    amrex::RealBox domain(AMREX_D_DECL(0.0,0.0,0.0), AMREX_D_DECL(1.0,1.0,1.0));
    AmrSim lbrx(nx, ny, nz, tau, tau, periodicity, domain);
    lbrx.SetInitialDensity(rho_0);
    lbrx.SetInitialVelocity(u_0);
    lbrx.InitFromScratch(0.0);

    printDensity(lbrx);
    printVelocity(lbrx);
    lbrx.Iterate(100);
    lbrx.CalcHydroVars(0);
    printDensity(lbrx);
    printVelocity(lbrx);
    lbrx.Iterate(100);
    lbrx.CalcHydroVars(0);
    printDensity(lbrx);
    printVelocity(lbrx);
  }
    amrex::Finalize();

  return 0;
}

void printDensity(const AmrSim& sim) {
  const std::array<int,NDIMS> DIMS = sim.GetDims();
  const int NUMEL = DIMS[0] * DIMS[1] * DIMS[2];
  const int LEVEL = 0;
  std::cout.precision(6);

  std::cout << "Time: " << sim.GetTime(LEVEL) << std::endl;
  std::cout << "Density:\n{ ";
  for (int k = 0; k < DIMS[2]; ++k) {
    for (int j = 0; j < DIMS[1]; ++j) {
      for (int i = 0; i < DIMS[0]; ++i) {
        if (k*DIMS[0]*DIMS[1] + j*DIMS[0] + i < NUMEL-1) {
          std::cout << sim.GetDensity(i, j, k, LEVEL) << ", ";
        } else {
          std::cout << sim.GetDensity(i, j, k, LEVEL) << " };";
        }
      }
    }
  }
  std::cout << std::endl;

  return;
}

void printVelocity(const AmrSim& sim) {
  const std::array<int,NDIMS> DIMS = sim.GetDims();
  const int NUMEL = DIMS[0] * DIMS[1] * DIMS[2] * NDIMS;
  std::cout.precision(6);
  const int LEVEL = 0;

  std::cout << "Time step: " << sim.GetTimeStep(LEVEL) << std::endl;
  std::cout << "Velocity:\n{ ";
  for (int k = 0; k < DIMS[2]; ++k) {
    for (int j = 0; j < DIMS[1]; ++j) {
      for (int i = 0; i < DIMS[0]; ++i) {
        for (int n = 0; n < NDIMS; ++n) {
          if (n*DIMS[0]*DIMS[1]*DIMS[2] + k*DIMS[0]*DIMS[1] + j*DIMS[0] + i < NUMEL-1) {
            std::cout << sim.GetVelocity(i, j, k, n, LEVEL) << ", ";
          } else {
            std::cout << sim.GetVelocity(i, j, k, n, LEVEL) << " };";
          }
        }
      }
    }
  }
  std::cout << std::endl;

  return;
}
