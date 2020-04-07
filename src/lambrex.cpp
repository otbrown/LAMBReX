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
