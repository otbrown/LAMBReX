#include "AmrTest.h"

void AmrTest::CallErrorEst(const int level, amrex::TagBoxArray& tba) {
  const double time = GetTime(level);
  const int ngrow = 1;

  ErrorEst(level, tba, time, ngrow);

  return;
}

void AmrTest::CallMakeNewLevelFromScratch(const amrex::BoxArray& ba,
const amrex::DistributionMapping& dm, const double time) {
  // MakeNewLevelFromScratch really only works at level 0,
  // happily this is also the only place it makes sense...
  const int level = 0;

  MakeNewLevelFromScratch(level, time, ba, dm);

  return;
}

void AmrTest::CallMakeNewLevelFromCoarse(const int level,
const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
  const double time = GetTime(level-1);

  MakeNewLevelFromCoarse(level, time, ba, dm);

  return;
}

void AmrTest::CallRemakeLevel(const int level, const double time,
const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
  RemakeLevel(level, time, ba, dm);
  return;
}

void AmrTest::CallClearLevel(const int level) {
  ClearLevel(level);
  return;
}
