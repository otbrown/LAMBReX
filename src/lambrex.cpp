#include <cstdio>
#include "AMReX.H"
#include "lambrex.h"

Simulation * lambrexInit(int argc, char ** argv, int nx, int ny, int nz, double tau_s, double tau_b, int (&periodicity)[NDIMS]) {
  printf("Initialising LAMBReX\n");

  amrex::Initialize(argc, argv, false);

  Simulation * lbrx = new Simulation(nx, ny, nz, tau_s, tau_b, periodicity);
  printf("A lamb was born!\n");

  return lbrx;
}

void lambrexFinalize(Simulation * lbrx) {
  printf("Finalising LAMBReX\n");

  // Simulation must be deleted before amrex::Finalize() is called,
  // as amrex::Finalize() invalidates some internal destructors
  delete lbrx;

  amrex::Finalize();

  return;
}
