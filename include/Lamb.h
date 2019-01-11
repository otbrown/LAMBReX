#ifndef LAMB_H
#define LAMB_H

#include "AMReX_Config.H"
#include "AMReX_Box.H"
#include "AMReX_RealBox.H"
#include "AMReX_IntVect.H"
#include "AMReX_IndexType.H"
#include "AMReX_Geometry.H"
#include "AMReX_FArrayBox.H"

// number of spatial dimensions and velocities are model dependent
// and cannot be changed (for now)
#define NDIMS 3
#define NMODES 15

class Lamb {
private:
  // physical dimensions of the outer domain
  const int NX;
  const int NY;
  const int NZ;
  const int NUMEL = NX * NY * NZ;
  const int COORD_SYS;
  const int PERIODICITY[NDIMS];

  // model parameters
  const double TAU_S; // shear
  const double TAU_B; // bulk
  const double CS2 = 1.0 / 3.0; // speed of sound squared
  static const double IDENTITY[NDIMS][NDIMS];
  static const double MODE_MATRIX[NMODES][NMODES];
  static const double MODE_MATRIX_INVERSE[NMODES][NMODES];

  // current time step
  int time_step = 0;

  //  AMReX domain specification
  amrex::Box idx_domain;
  amrex::RealBox phys_domain;
  amrex::Geometry geometry;

  // hydrodynamic variables (output arrays)
  amrex::FArrayBox density;
  amrex::FArrayBox velocity;
  amrex::FArrayBox force;

  // distribution function (work array)
  amrex::FArrayBox dist_fn;

  // member functions
  void collide();
  void propagate();

public:
  Lamb(int, int, int, double, double, int (&)[NDIMS]);

  int getTimeStep() { return time_step; }
  void setDensity(double);
  void setDensity(double *);
  void setVelocity(double);
  void setVelocity(double *);
  void calcEquilibriumDist();

  int iterate(int);

  void printDensity();
};

#endif
