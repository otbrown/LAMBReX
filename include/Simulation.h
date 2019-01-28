#ifndef SIMULATION_H
#define SIMULATION_H

#include "AMReX_Config.H"
#include "AMReX_Box.H"
#include "AMReX_RealBox.H"
#include "AMReX_IntVect.H"
#include "AMReX_IndexType.H"
#include "AMReX_Geometry.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_BoxArray.H"
#include "AMReX_DistributionMapping.H"
#include "AMReX_MultiFab.H"

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
  constexpr static double CS2 = 1.0 / 3.0; // speed of sound squared
  constexpr static int HALO_DEPTH = 1;
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
  amrex::BoxArray ba_domain;
  amrex::DistributionMapping dm;

  // hydrodynamic variables (output arrays)
  amrex::FArrayBox density;
  amrex::FArrayBox velocity;
  amrex::FArrayBox force;

  // distribution function (work array)
  amrex::MultiFab dist_fn;

  // member functions
  int lindex(const int i, const int j, const int k, const int n) const {
    return (n * (NX+2*HALO_DEPTH) * (NY+2*HALO_DEPTH) * (NZ+2*HALO_DEPTH)
            + k * (NX+2*HALO_DEPTH) * (NY+2*HALO_DEPTH) + j * (NX+2*HALO_DEPTH)
            + i);
  }
  void swapElements(double * const, const int, const int);
  void updateBoundaries();
  void collide();
  void propagate();

public:
  Simulation(int const, int const, int const, double const, double const,
    int (&)[NDIMS]);

  int getTimeStep() const { return time_step; }
  std::array<int,NDIMS> getDims() const {return std::array<int,NDIMS>{NX,NY,NZ};}
  double getDensity(const int, const int, const int) const;
  void setDensity(double const);
  void setDensity(double const * const);
  void setVelocity(double const);
  void setVelocity(double const * const);
  void calcEquilibriumDist();
  void calcHydroVars();

  int iterate(int const);

  void printDensity() const;
};

#endif
