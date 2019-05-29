#ifndef AMRSIM_H
#define AMRSIM_H

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"
#include "AMReX_BCRec.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_Interpolater.H"

// number of spatial dimensions and velocities are model dependent
// and cannot be changed (for now)
#define NDIMS 3
#define NMODES 15

class AmrSim : public amrex::AmrCore {
protected:
  // physical dimensions of the outer domain
  const int NX;
  const int NY;
  const int NZ;
  const int NUMEL;
  const int COORD_SYS;
  std::array<int,NDIMS> PERIODICITY;

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

  // boundary conditions (in AMReX format)
  // at the moment we only allow periodic boundary conditions
  amrex::Vector<amrex::BCRec> f_bndry;
  const amrex::BndryFuncArray bfunc;

  // address of interpolator to use between coarse and fine grids
  // these are constructed as globals by AMReX...
  // cell_cons_interp : cell conservative linear interpolation
  amrex::Interpolater * mapper = &amrex::cell_cons_interp;

  // current time steps
  std::vector<double> sim_time;
  std::vector<double> dt;
  std::vector<int> time_step;

  // hydrodynamic variables (output arrays)
  std::vector<amrex::MultiFab> density;
  std::vector<amrex::MultiFab> velocity;

  // distribution function (work array)
  std::vector<amrex::MultiFab> dist_fn;

  // member functions
  int FLindex(int const i, int const j, int const k, int const n,
              amrex::IntVect const dims) const {
    return (n * dims[0] * dims[1] * dims[2] + k * dims[0] * dims[1]
            + j * dims[0] + i);
  }
  int CLindex(int const i, int const j, int const k, int const n,
              amrex::IntVect const dims, int const n_comps) const {
    return (i * dims[1] * dims[2] * n_comps + j * dims[2] * n_comps
            + k * n_comps + n);
  }
  void SwapElements(double * const, int const, int const);
  void UpdateBoundaries(int const);
  void Collide(int const);
  void Propagate(int const);
  void InitDensity(int const);
  void InitVelocity(int const);
  void ComputeDt(int const);
  void IterateLevel(int const);
  static void DistFnFillShim(double *, const int *, const int *, const int *, const int *,
    const double *, const double *, const double *, const int *);
  void DistFnFillPatch(int const, amrex::MultiFab&);
  void DistFnFillFromCoarse(int const, amrex::MultiFab&);
  bool TagCell(const amrex::FArrayBox&, const amrex::IntVect&);

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
  AmrSim(double const, double const);
  double GetTime(int const level) const { return sim_time.at(level); }
  int GetTimeStep(int const level) const { return time_step.at(level); }
  std::array<int,NDIMS> GetDims() const {
    return std::array<int,NDIMS>{NX,NY,NZ}; }
  void SetInitialDensity(double const);
  void SetInitialDensity(std::vector<double> const);
  void SetInitialVelocity(double const);
  void SetInitialVelocity(std::vector<double> const);
  double GetDensity(int const, int const, int const, int const) const;
  double GetVelocity(int const, int const, int const, int const, int const)
    const;
  void CalcEquilibriumDist(int const);
  void CalcHydroVars(int const);
  void Iterate(int const);
};

#endif
