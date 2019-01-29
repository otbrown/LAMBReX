#include <iostream>
#include <memory>
#include "lambrex.h"

void report(const Simulation&);

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
    lbrx.setDensity(rho);
    std::cout << "Density initialised." << std::endl;
    lbrx.setVelocity(0.0);
    std::cout << "Velocity initialised." << std::endl;

    lbrx.calcEquilibriumDist();
    std::cout << "Equilibrium distribution calculated." << std::endl;

    report(lbrx);
    lbrx.iterate(100);
    lbrx.calcHydroVars();
    report(lbrx);
    lbrx.iterate(100);
    lbrx.calcHydroVars();
    report(lbrx);
  }
  amrex::Finalize();

  return 0;
}

void report(const Simulation& sim) {
  const std::array<int,NDIMS> DIMS = sim.getDims();
  std::cout.precision(6);

  std::cout << "Time step: " << sim.getTimeStep() << std::endl;
  std::cout << "Density:\n[ ";
  for (int k = 0; k < DIMS[2]; ++k) {
    std::cout << sim.getDensity(DIMS[0]/2, DIMS[1]/2, k) << " ";
  }
  std::cout << "]" << std::endl;

  return;
}
