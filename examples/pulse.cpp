#include <iostream>
#include "lambrex.h"

void printDensity(const Simulation&);
void printVelocity(const Simulation&);

int main (int argc, char * argv[])
{
  const int nx = 10;
  const int ny = 10;
  const int nz = 50;
  const double tau = 0.5;
  double amplitude = 0.01;
  int periodicity[3] = {1, 1, 1};
  const int numel = nx * ny * nz;
  int i, j, k;
  double z_mean;

  // create initial density array
  double rho[numel] = {};
  for (i = 0; i < numel; ++i) rho[i] = 1.0;

  // set enhanced density
  for (i = 0; i < nx; ++i) {
    for (j = 0; j < ny; ++j) {
      rho[i*ny*nz + j*nz + nz/2 - 1] += amplitude;
    }
  }

  // renormalise
  z_mean = 0.0;
  for (k = 0; k < nz; ++k) z_mean += rho[((numel + ny * nz) / 2) + k];
  z_mean /= nz;
  for (i = 0; i < numel; ++i) rho[i] /= z_mean;

  amrex::Initialize(argc, argv, false);
  {
    Simulation lbrx(nx, ny, nz, tau, tau, periodicity);
    // provide LAMBReX with initial density and velocity
    // density is copied so safe to free
    for (i = 0; i < nx; ++i) {
      for (j = 0; j < ny; ++j) {
        for (k = 0; k < nz; ++k) {
          lbrx.setDensity(i, j, k, rho[i*ny*nz + j*nz + k]);
        }
      }
    }
    std::cout << "Density initialised." << std::endl;
    lbrx.setVelocity(0.0);
    std::cout << "Velocity initialised." << std::endl;

    lbrx.calcEquilibriumDist();
    std::cout << "Equilibrium distribution calculated." << std::endl;

    printDensity(lbrx);
    printVelocity(lbrx);
    lbrx.iterate(100);
    lbrx.calcHydroVars();
    printDensity(lbrx);
    printVelocity(lbrx);
    lbrx.iterate(100);
    lbrx.calcHydroVars();
    printDensity(lbrx);
    printVelocity(lbrx);
  }
  amrex::Finalize();

  return 0;
}

void printDensity(const Simulation& sim) {
  const std::array<int,NDIMS> DIMS = sim.getDims();
  const int NUMEL = DIMS[0] * DIMS[1] * DIMS[2];
  std::cout.precision(6);

  std::cout << "Time step: " << sim.getTimeStep() << std::endl;
  std::cout << "Density:\n{ ";
  for (int k = 0; k < DIMS[2]; ++k) {
    for (int j = 0; j < DIMS[1]; ++j) {
      for (int i = 0; i < DIMS[0]; ++i) {
        if (k*DIMS[0]*DIMS[1] + j*DIMS[0] + i < NUMEL-1) {
          std::cout << sim.getDensity(i, j, k) << ", ";
        } else {
          std::cout << sim.getDensity(i, j, k) << " };";
        }
      }
    }
  }
  std::cout << std::endl;

  return;
}

void printVelocity(const Simulation& sim) {
  const std::array<int,NDIMS> DIMS = sim.getDims();
  const int NUMEL = DIMS[0] * DIMS[1] * DIMS[2] * NDIMS;
  std::cout.precision(6);

  std::cout << "Time step: " << sim.getTimeStep() << std::endl;
  std::cout << "Velocity:\n{ ";
  for (int k = 0; k < DIMS[2]; ++k) {
    for (int j = 0; j < DIMS[1]; ++j) {
      for (int i = 0; i < DIMS[0]; ++i) {
        for (int n = 0; n < NDIMS; ++n) {
          if (n*DIMS[0]*DIMS[1]*DIMS[2] + k*DIMS[0]*DIMS[1] + j*DIMS[0] + i < NUMEL-1) {
            std::cout << sim.getVelocity(i, j, k, n) << ", ";
          } else {
            std::cout << sim.getVelocity(i, j, k, n) << " };";
          }
        }
      }
    }
  }
  std::cout << std::endl;

  return;
}
