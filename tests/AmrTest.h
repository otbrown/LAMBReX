#ifndef AMRTEST_H
#define AMRTEST_H

#include "AmrSim.h"

class AmrTest : public AmrSim {
public:
  using AmrSim::AmrSim;
  using AmrSim::GetDensity;
  using AmrSim::GetVelocity;
  using AmrSim::GetTimeStep;

  bool DensityEmpty(const int level) const { return density.at(level).empty(); }
  bool VelocityEmpty(const int level) const { return velocity.at(level).empty(); }
  bool DistFnEmpty(const int level) const { return dist_fn.at(level).empty(); }

  std::vector<amrex::MultiFab>& GetDensity() { return density; };
  std::vector<amrex::MultiFab>& GetVelocity() { return velocity; };
  std::vector<amrex::MultiFab>& GetDistFn() { return dist_fn; };
  std::vector<double>& GetSimTime() { return sim_time; }
  std::vector<double>& GetDt() { return dt; }
  std::vector<int>& GetTimeStep() { return time_step; }
  std::vector<double>& GetTauS() { return tau_s; }
  std::vector<double>& GetTauB() { return tau_b; }
  std::vector<double>& GetMass() { return mass; }

  void CallErrorEst(int const, amrex::TagBoxArray&);
  void CallMakeNewLevelFromScratch(const int, const amrex::BoxArray&,
                                   const amrex::DistributionMapping&,
                                   const double);
  void CallMakeNewLevelFromCoarse(const int, const amrex::BoxArray&,
                                  const amrex::DistributionMapping&);
  void CallRemakeLevel(const int, const double, const amrex::BoxArray&,
                       const amrex::DistributionMapping&);
  void CallClearLevel(const int);
};

#endif
