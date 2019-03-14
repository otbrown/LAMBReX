#include "lambrex.h"

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

  amrex::Initialize(argc, argv, false);
  {
    amrex::RealBox domain(AMREX_D_DECL(0.0,0.0,0.0), AMREX_D_DECL(1.0,1.0,1.0));
    AmrSim lbrx(nx, ny, nz, tau, tau, periodicity, domain);
  }
  amrex::Finalize();

  return 0;
}
