#ifndef LAMB_H
#define LAMB_H

class Lamb {
private:
  // physical dimensions of the outer domain
  const int _NDIMS = 3;
  const int _NX;
  const int _NY;
  const int _NZ;

  // relaxation times
  const double _TAU_S; // shear
  const double _TAU_B; // bulk

  // speed of sound squared
  const double _CS2 = 1.0 / 3.0;

  // current time step
  int time_step = 0;

public:
  Lamb(int, int, int, double, double);
};

#endif
