#ifndef AMRSIM_H
#define AMRSIM_H

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"

// number of spatial dimensions and velocities are model dependent
// and cannot be changed (for now)
#define NDIMS 3
#define NMODES 15

class AmrSim : public amrex::AmrCore {
private:
  // physical dimensions of the outer domain
  const int NX;
  const int NY;
  const int NZ;
  const int NUMEL;
  const int COORD_SYS;
  int PERIODICITY[NDIMS];

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
  std::vector<double> initial_density;
  std::vector<double> initial_velocity;

  // current time step
  std::vector<double> time_step;

  //  AMReX domain specification
  amrex::Box idx_domain;
  amrex::RealBox phys_domain;
  amrex::BoxArray ba_domain;
  amrex::DistributionMapping dm;

  // hydrodynamic variables (output arrays)
  std::vector<amrex::MultiFab> density;
  std::vector<amrex::MultiFab> velocity;

  // distribution function (work array)
  std::vector<amrex::MultiFab> dist_fn;

  // member functions
  int Lindex(const int i, const int j, const int k, const int n,
             const amrex::IntVect dims) const {
    return (n * dims[0] * dims[1] * dims[2] + k * dims[0] * dims[1]
            + j * dims[0] + i);
  }
  void SwapElements(double * const, const int, const int);
  void UpdateBoundaries(int);
  void Collide(int);
  void Propagate(int);

  // AMRCore pure virtual functions
  void ErrorEst(int, amrex::TagBoxArray&, double, int) override;
  void MakeNewLevelFromScratch(int, double, const amrex::BoxArray&,
    const amrex::DistributionMapping&) override;
  void MakeNewLevelFromCoarse(int, double, const amrex::BoxArray&,
    const amrex::DistributionMapping&) override;
  void RemakeLevel(int, double, const amrex::BoxArray&,
    const amrex::DistributionMapping&) override;
  void ClearLevel(int) override;

public:
  AmrSim(int const, int const, int const, double const, double const,
    int (&)[NDIMS], amrex::RealBox&);
  double GetTimeStep(int level) const { return time_step.at(level); }
  std::array<int,NDIMS> GetDims() const {
    return std::array<int,NDIMS>{NX,NY,NZ}; }
  void SetInitialDensity(const double);
  void SetInitialDensity(const std::vector<double>);
  void SetInitialVelocity(const double);
  void SetInitialVelocity(const std::vector<double>);
  void InitDistFunc();
  double GetDensity(int const, int const, int const, int const) const;
  double GetVelocity(int const, int const, int const, int const, int const)
    const;
  void CalcEquilibriumDist();
  void CalcHydroVars();

  int Iterate(int const);

  void PrintDensity() const;
};

#endif
