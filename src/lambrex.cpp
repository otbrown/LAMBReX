#include "lambrex.h"
#include "AMReX_ParmParse.H"

void lambrexInit() {
  int argc = 0;
  char ** argv;
  amrex::Initialize(argc, argv, false);
  return;
}

void lambrexFinalise() {
  amrex::Finalize();
  return;
}
