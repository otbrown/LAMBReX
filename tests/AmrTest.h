#ifndef AMRTEST_H
#define AMRTEST_H

#include "AmrSim.h"

class AmrTest : public AmrSim {
  using AmrSim::AmrSim;

public:
  bool DensityEmpty(const int level) { return density.at(level).empty(); }
  bool VelocityEmpty(const int level) { return velocity.at(level).empty(); }
  bool DistFnEmpty(const int level) { return dist_fn.at(level).empty(); }
  void CallErrorEst(int const, amrex::TagBoxArray&);
  void CallMakeNewLevelFromScratch(const amrex::BoxArray&,
                                   const amrex::DistributionMapping&,
                                   const double);
  void CallMakeNewLevelFromCoarse(const int, const amrex::BoxArray&,
                                  const amrex::DistributionMapping&);
  void CallRemakeLevel(const int, const double, const amrex::BoxArray&,
                       const amrex::DistributionMapping&);
  void CallClearLevel(const int);
};

#endif
