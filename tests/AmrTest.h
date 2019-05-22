#ifndef AMRTEST_H
#define AMRTEST_H

#include "AmrSim.h"

class AmrTest : public AmrSim {
public:
  void CallErrorEst(int const, amrex::TagBoxArray&);
  void CallMakeNewLevelFromScratch(const amrex::BoxArray&,
                                   const amrex::DistributionMapping&);
  void CallMakeNewLevelFromCoarse(const int level, const amrex::BoxArray&,
                                  const amrex::DistributionMapping&);
  void CallRemakeLevel(const int, const double, const amrex::BoxArray&,
                       const amrex::DistributionMapping&);
  void CallClearLevel(const int);
};

#endif
