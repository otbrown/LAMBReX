#ifndef LAMB_H
#define LAMB_H

#include "AMReX_Config.H"
#include "AMReX_Box.H"
#include "AMReX_RealBox.H"
#include "AMReX_Geometry.H"

class Lamb {
private:
  // physical dimensions of the outer domain
  const int _NDIMS = AMREX_SPACEDIM;
  const int _NX;
  const int _NY;
  const int _NZ;
  const int _NUMEL = _NX * _NY * _NZ;
  const int COORD_SYS = 0;
  int _periodicity[AMREX_SPACEDIM];

  // relaxation times
  const double _TAU_S; // shear
  const double _TAU_B; // bulk

  // speed of sound squared
  const double CS2 = 1.0 / 3.0;

  // current time step
  int time_step = 0;

  //  AMReX domain specification
  amrex::Box idx_domain;
  amrex::RealBox phys_domain;
  amrex::Geometry geometry;

  // hydrodynamic variables (output arrays)
  double * density;
  double * velocity;
  double * force;

  void buildGeo();

public:
  Lamb(int, int, int, double, double, int *);
  ~Lamb();
};

#endif
