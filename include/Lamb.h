#ifndef LAMB_H
#define LAMB_H

#include "AMReX_Config.H"
#include "AMReX_Box.H"
#include "AMReX_RealBox.H"
#include "AMReX_Geometry.H"
#include "AMReX_FArrayBox.H"

class Lamb {
private:
  // physical dimensions of the outer domain
  const int NDIMS = AMREX_SPACEDIM;
  const int NX;
  const int NY;
  const int NZ;
  const int NUMEL = NX * NY * NZ;
  const int COORD_SYS = 0;
  int PERIODICITY[AMREX_SPACEDIM];

  // model parameters
  const int NUM_VELOCITIES = 15;
  const double TAU_S; // shear
  const double TAU_B; // bulk
  // speed of sound squared
  const double CS2 = 1.0 / 3.0;

  // current time step
  int time_step = 0;

  //  AMReX domain specification
  amrex::Box idx_domain;
  amrex::RealBox phys_domain;
  amrex::Geometry geometry;

  // hydrodynamic variables (output arrays)
  amrex::FArrayBox * density;
  amrex::FArrayBox * velocity;
  amrex::FArrayBox * force;

  // distribution function (work array)
  amrex::FArrayBox * dist_fn;

  // member functions
  void buildGeometry();

public:
  Lamb(int, int, int, double, double, int *);
  ~Lamb();

  void setDensity(double);
  void setDensity(double *);
  void setVelocity(double);
  void setVelocity(double *);
  void calcEquilibriumDist();
};

#endif
