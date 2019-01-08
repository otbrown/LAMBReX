#include <cstdio>
#include "lambrex.h"

int main (int argc, char * argv[])
{
  int nx = 10;
  int ny = 10;
  int nz = 50;
  double tau = 0.5;
  int periodicity[3] = {1, 1, 1};
  int numel = nx * ny * nz;
  double amplitude = 0.01;
  int i, j, k;
  double z_mean;

  // create initial density array
  double * rho = new double[numel];
  for (i = 0; i < numel; ++i) rho[i] = 1.0;
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

  Lamb * lbrx = lambrexInit(argc, argv, nx, ny, nz, tau, tau, periodicity);

  // provide LAMBReX with initial density and velocity
  // density is memcpy'd so safe to free
  lbrx->setDensity(rho);
  printf("Density initialised.\n");
  lbrx->setVelocity(0.0);
  printf("Velocity initialised.\n");

  delete[] rho;

  lambrexFinalize(lbrx);

  return 0;
}
