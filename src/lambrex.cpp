#include <cstdio>
#include "lambrex.h"

Lamb * lambrexInit(int nx, int ny, int nz, double tau_s, double tau_b) {
  printf("Initialising LAMBReX\n");

  Lamb * lbrx = new Lamb(nx, ny, nz, tau_s, tau_b);
  printf("A lamb was born!\n");

  return lbrx;
}

void lambrexFinalize(Lamb * lbrx) {
  printf("Finalising LAMBReX\n");

  delete lbrx;
  printf("Lamb deleted.\n");

  return;
}
