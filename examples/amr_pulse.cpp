#include <iostream>
#include "lambrex.h"

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
    std::cout << "Initialised lbrx" << std::endl;
    lbrx.SetInitialDensity(rho_0);
    std::cout << "Set initial density" << std::endl;
    lbrx.SetInitialVelocity(u_0);
    std::cout << "Set initial velocity" << std::endl;
    lbrx.InitFromScratch(0.0);
    std::cout << "Initialised lbrx MultiFabs" << std::endl;
  }
  amrex::Finalize();

  return 0;
}
