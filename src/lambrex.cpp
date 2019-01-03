#include <cstdio>
#include "AMReX.H"
#include "lambrex.h"

Lamb * lambrexInit(int argc, char ** argv, int nx, int ny, int nz, double tau_s, double tau_b, int * periodicity) {
  printf("Initialising LAMBReX\n");

  amrex::Initialize(argc, argv, false);

  Lamb * lbrx = new Lamb(nx, ny, nz, tau_s, tau_b, periodicity);
  printf("A lamb was born!\n");

  return lbrx;
}

void lambrexFinalize(Lamb * lbrx) {
  printf("Finalising LAMBReX\n");

  amrex::Finalize();

  delete lbrx;
  printf("Lamb deleted.\n");

  return;
}
