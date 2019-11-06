// -*- mode: c++; -*-
#ifndef AMRSIM_H
#define AMRSIM_H

#include "AMReX_AmrCore.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFab.H"
#include "AMReX_BCRec.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_Interpolater.H"

#include "velocity_set.h"
#include "component.h"
#include "multilevel.h"
#include "d3q15_bgk.h"

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
  constexpr static double NL_DENSITY = -1.0;
  constexpr static double NL_VELOCITY = -3E8;
  static const double DELTA[NDIMS][NDIMS];
  static const double MODE_MATRIX[NMODES][NMODES];
  static const double MODE_MATRIX_INVERSE[NMODES][NMODES];
  std::vector<double> initial_density;
  std::vector<double> initial_velocity;
  std::vector<amrex::BoxArray> static_tags;

  // boundary conditions (in AMReX format)
  // at the moment we only allow periodic boundary conditions
  amrex::Vector<amrex::BCRec> f_bndry;
  const amrex::BndryFuncArray bfunc;

  // address of interpolator to use between coarse and fine grids
  // these are constructed as globals by AMReX...
  // cell_cons_interp : cell conservative linear interpolation
  amrex::Interpolater * mapper = &amrex::cell_cons_interp;

  std::vector<double> tau_s;
  std::vector<double> tau_b;

  // mass constant per level, may be needed for calculating outputs
  std::vector<double> mass;

  // hydrodynamic variables (output arrays)
  std::vector<amrex::MultiFab> velocity;

  // distribution function (work array)
  std::vector<SimLevelData> levels;

  // masks to keep track of which parts of the coarse level are covered by a
  // fine grid. The iMultiFab contains COARSE_VAL if not, FINE_VAL if so. The
  // convention from AMReX itself is that COARSE_VAL=0, FINE_VAL=1.
  const int COARSE_VAL = 0;
  const int FINE_VAL = 1;
  std::vector<amrex::iMultiFab> fine_masks;

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
  void Collide(const amrex::MultiFab&, amrex::MultiFab&, const double, const double);
  void Stream(int const);
  void CollideLevel(int const);
  void CollideAndStream(int const LEVEL) {
    CollideLevel(LEVEL);
    Stream(LEVEL);
    levels.at(LEVEL).UpdateNow();
    return;
  }
  void InitDensity(int const);
  void InitVelocity(int const);
  void ComputeDt(int const);
  void IterateLevel(int const);
  void SubCycle(int const, int const);
  static void DistFnFillShim(double *, const int *, const int *, const int *,
    const int *, const double *, const double *, const double *, const int *);
  void DistFnFillPatch(int const, amrex::MultiFab&);
  void DistFnFillFromCoarse(int const, amrex::MultiFab&);
  bool TagCell(int const, const amrex::IntVect&);
  void MakeFineMask(int const);

  // Rohde Steps
  void RohdeCycle(int const);
  void Explode(int const);
  void SumFromFine(int const);
  void StreamInterior(const amrex::MultiFab&, amrex::MultiFab&);
  void CompleteTimeStep(int const);

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
  inline double GetTime(int const level) const {
    return levels.at(level).time.current;
  }
  int GetTimeStep(int const level) const {
    return levels.at(level).time.step;
  }
  std::array<int,NDIMS> GetDims() const {
    return std::array<int,NDIMS>{NX,NY,NZ}; }
  bool OnProcessDensity(double const rho) const { return rho != NL_DENSITY; }
  bool OnProcessVelocity(double const u) const { return u != NL_VELOCITY; }
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
  void SetStaticRefinement(int const, const std::array<int, NDIMS>&,
    const std::array<int, NDIMS>&);
  void UnsetStaticRefinement(int const);
  std::pair<std::array<int,NDIMS>, std::array<int,NDIMS>> GetExtent(int const)
  const;
};

#endif
