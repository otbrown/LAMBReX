#include "Lamb.h"

Lamb::Lamb(int nx, int ny, int nz, double tau_s, double tau_b) : _NX(nx), _NY(ny), _NZ(nz), _TAU_S(tau_s), _TAU_B(tau_b) {
  density = new double[_NUMEL];
  velocity = new double[_NUMEL];
  force = new double[_NUMEL];

  return;
};

Lamb::~Lamb() {
  delete[] density;
  delete[] velocity;
  delete[] force;

  return;
}
