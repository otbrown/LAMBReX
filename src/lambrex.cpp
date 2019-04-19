#include "lambrex.h"
#include "AMReX_ParmParse.H"

void lambrexInit(const std::array<int,NDIMS>& periodicity) {
  int argc = 0;
  char ** argv;

  {
    amrex::ParmParse pp;
    pp.add("title", "AmrSim");
  }
  {
    amrex::ParmParse pp("amr");
    pp.add("v", 1);
  }
  // periodicity is set here as, despite appearances, it is *not* configurable
  // later on, at least for level 0
  {
    amrex::ParmParse pp("geometry");
    pp.addarr("is_periodic",
      std::vector<int>(periodicity.begin(), periodicity.end()));
    pp.add("coord_sys", 0); // 0 -> Cartesian coordinates
    // determine "physical" problem domain size
    pp.addarr("prob_lo", std::vector<double>({0.0, 0.0, 0.0}));
    pp.addarr("prob_hi", std::vector<double>({1.0, 1.0, 1.0}));
  }

  amrex::Initialize(argc, argv, true);

  return;
}

void lambrexFinalise(void) {
  amrex::Finalize();
  return;
}

void lambrexSetAmr(const int nx, const int ny, const int nz,
const int max_ref_level) {
  const int blocking_factor = binGCD(nx, binGCD(ny, nz));
  std::cout << "Blocking factor set to: " << blocking_factor << std::endl;

  {
    amrex::ParmParse pp("amr");
    // n_cell determines no. of grid points in every direction
    pp.addarr("n_cell", std::vector<int>({nx, ny, nz}));
    pp.add("max_level", max_ref_level);
    pp.add("blocking_factor", blocking_factor);
    pp.add("max_grid_size", 2*blocking_factor);
  }

  return;
}

int binGCD(int a, int b) {
  // greatest common divisor using Stein's Method
  // std::gcd is C++17...
  if (a < 0 || b < 0) amrex::Abort("Negative domain sizes are NOT allowed!");

  // simple cases
  if (a == b) return a;
  if (!a) return b;
  if (!b) return a;

  if (!(a%2)) {
    if (!(b%2)) return 2 * binGCD(a/2, b/2);  // a and b even
    else return binGCD(a/2, b); // a even, b odd
  } else {
    if (!(b%2)) return binGCD(a, b/2);  // a odd, b even
  }

  if (a >= b) return binGCD((a-b)/2, b);  // a odd, b odd, a >= b
  else return binGCD((b-a)/2, a); // a odd, b odd, a < b
}
