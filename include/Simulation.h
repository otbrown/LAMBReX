#ifndef SIMULATION_H
#define SIMULATION_H

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

class Simulation {
private:
  // physical dimensions of the outer domain
  const int NX;
  const int NY;
  const int NZ;
  const int NUMEL;
  const int COORD_SYS;
  const int PERIODICITY[NDIMS];

  // model parameters
  const double CS2 = 1.0 / 3.0; // speed of sound squared
  const double TAU_S; // shear
  const double TAU_B; // bulk
  const double OMEGA_S;
  const double OMEGA_B;
  static const double DELTA[NDIMS][NDIMS];
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
  Simulation(int const, int const, int const, double const, double const,
    int (&)[NDIMS]);

  int getTimeStep() { return time_step; }
  void setDensity(double const);
  void setDensity(double const * const);
  void setVelocity(double const);
  void setVelocity(double const * const);
  void calcEquilibriumDist();

  int iterate(int const);

  void printDensity();
};

#endif
