#include <cstdio>
#include <memory>
#include "lambrex.h"

int main (int argc, char * argv[])
{
  const int nx = 10;
  const int ny = 10;
  const int nz = 50;
  const double tau = 0.5;
  int periodicity[3] = {1, 1, 1};
  const int numel = nx * ny * nz;
  double amplitude = 0.01;
  int i, j, k;
  double z_mean;

  // create initial density array
  double rho[numel] = {1.0};

  // set enhanced density
  for (i = 0; i < nx; ++i) {
    for (j = 0; j < ny; ++j) {
      rho[i*ny*nz + j*nz + nz/2] += amplitude;
    }
  }
  // renormalise
  z_mean = 0.0;
  for (k = 0; k < nz; ++k) z_mean += rho[((numel + ny * nz) / 2) + k];
  z_mean /= nz;
  for (i = 0; i < numel; ++i) rho[i] /= z_mean;

  amrex::Initialize(argc, argv, false);
  {
    const std::unique_ptr<Simulation> lbrx =
      std::make_unique<Simulation>(nx, ny, nz, tau, tau, periodicity);

    // provide LAMBReX with initial density and velocity
    // density is copied so safe to free
    lbrx->setDensity(rho);
    printf("Density initialised.\n");
    lbrx->setVelocity(0.0);
    printf("Velocity initialised.\n");

    lbrx->calcEquilibriumDist();
    printf("Equilibrium distribution calculated.\n");
  }
  amrex::Finalize();

  return 0;
}
