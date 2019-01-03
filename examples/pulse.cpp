#include <cstdio>
#include "lambrex.h"

int main (int argc, char * argv[])
{
  int nx = 10;
  int ny = 10;
  int nz = 50;
  double tau = 0.5;
  int periodicity[3] = {1, 1, 1};

  Lamb * lbrx = lambrexInit(argc, argv, nx, ny, nz, tau, tau, periodicity);

  printf("Hello world.\n");

  lambrexFinalize(lbrx);

  return 0;
}
