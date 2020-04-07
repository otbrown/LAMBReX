#ifndef AMRTEST_H
#define AMRTEST_H

#include "AmrSim.h"

class AmrTest : public AmrSim {
public:
  using AmrSim::AmrSim;
  using AmrSim::GetDensity;
  using AmrSim::GetVelocity;
  using AmrSim::GetTimeStep;

  bool DensityEmpty(const int level) const { return levels[level].now.get<Density>().empty(); }
  bool VelocityEmpty(const int level) const { return velocity.at(level).empty(); }
  bool DistFnEmpty(const int level) const { return levels[level].now.get<DistFn>().empty(); }

  const int CoarseVal() { return COARSE_VAL; }
  const int FineVal() { return FINE_VAL; }
  std::vector<amrex::MultiFab>& GetVelocity() { return velocity; };
  inline auto& GetLevels() {
    return levels;
  }
  inline auto& GetLevel(int l) {
    return levels[l];
  }

  inline auto& GetSimTime(int l) {
    return levels[l].time.current;
  }
  inline auto& GetDt(int l) {
    return levels[l].time.delta;
  }
  inline auto& GetTimeStep(int l) {
    return levels[l].time.step;
  }
  std::vector<double>& GetTauS() { return tau_s; }
  std::vector<double>& GetTauB() { return tau_b; }
  std::vector<double>& GetMass() { return mass; }
  amrex::iMultiFab& GetFineMask(int const LEVEL) {return fine_masks.at(LEVEL);}

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
