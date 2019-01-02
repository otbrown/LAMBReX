#ifndef LAMB_H
#define LAMB_H

#include "AMReX_Config.H"
#include "AMReX_Box.H"

class Lamb {
private:
  // physical dimensions of the outer domain
  const int _NDIMS = 3;
  const int _NX;
  const int _NY;
  const int _NZ;
  const int _NUMEL = _NX * _NY * _NZ;

  // relaxation times
  const double _TAU_S; // shear
  const double _TAU_B; // bulk

  // speed of sound squared
  const double _CS2 = 1.0 / 3.0;

  // current time step
  int time_step = 0;

  // distribution function
  amrex::Box dist_f;

  // hydrodynamic variables (output arrays)
  double * density;
  double * velocity;
  double * force;

public:
  Lamb(int, int, int, double, double);
  ~Lamb();
};

#endif
