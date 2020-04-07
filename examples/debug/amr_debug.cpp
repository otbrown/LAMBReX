#include <iostream>
#include "lambrex.h"
#include "AmrTest.h"

void printDensity(AmrTest&, const int);

int main (void) 
{
  // simulation set-up, cubic grid
  // refinement covering interior
  const int N = 8;
  const int NUMEL = N * N * N;
  std::array<int,3> periodicity = {1, 1, 1};
  const int MAX_LEVEL = 1;
  // REF : "REFINEMENT"
  // refinement grid defined by upper and lower corner
  const int REF_OFFSET = 2;
  const std::array<int,3> REF_LO{REF_OFFSET, REF_OFFSET, REF_OFFSET};
  const std::array<int,3> REF_HI{
    N-(REF_OFFSET+1), N-(REF_OFFSET+1), N-(REF_OFFSET+1)
    };
  
  // initial density {0, 1, 2, 3, ... N^3-1}
  double rho_i = 0.0;
  std::vector<double> rho_0(NUMEL);
  for (auto& density : rho_0) {
    density = rho_i;
    rho_i += 1.0;
  }

  // initial velocity 
  double u_0 = 0.0;

  // relaxation parameter
  const double TAU = 0.5;

  // initialise AmrSim
  lambrexInit(periodicity);
  lambrexSetAmr(N, N, N, MAX_LEVEL);
  {
    // AmrTest is a derived class which has public functions exposing internals
    AmrTest sim(TAU, TAU);
    sim.SetInitialDensity(rho_0);
    sim.SetInitialVelocity(u_0);
    sim.InitFromScratch(0.0);
    sim.SetStaticRefinement(0, REF_LO, REF_HI);

    sim.CalcHydroVars(0);
    printDensity(sim, 0);

    sim.Iterate(1);
    sim.CalcHydroVars(0);
    printDensity(sim, 0);
   
  }
  lambrexFinalise();

  return 0;
}

void printDensity(AmrTest& sim, const int LEVEL) {
  const auto& RHO_mf = sim.GetLevel(LEVEL).now.get<Density>();
  amrex::IntVect lo(0);
  amrex::IntVect hi(0);
  amrex::IntVect pos(0);

  std::cout << "Density on level " << LEVEL << ": " << "\n\n";
  int box_count = 0;
  for (amrex::MFIter mfi(RHO_mf); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const amrex::FArrayBox& rho = RHO_mf[mfi];

    lo = box.smallEnd();
    hi = box.bigEnd();

    std::cout << "Box " << box_count << ":\n\n";
    for (int k = lo[2]; k<= hi[2]; ++k) {
      pos.setVal(2, k);
      std::cout << "k = " << k << "\n\n";
      for (int i = lo[0]; i <= hi[0]; ++i) {
        pos.setVal(0, i);
        for(int j = lo[1]; j <= hi[1]; ++j) {
          pos.setVal(1, j);
          std::cout << rho(pos) << " ";
        }
        std::cout << "\n";
      }
      std::cout << "\n"; 
    }

    ++box_count;
  }

  return;
}