// Do some basic compile-time tests of the metaprogramming.
//
// If this compiles, then the tests pass.
//
#include "d3q15_bgk.h"
#include "range.hpp"

// Velocity set
static_assert(D3Q15::ND == 3);
static_assert(D3Q15::NV == 15);
static_assert(D3Q15::HALO == 2);

constexpr bool check_vels() {
  auto ans = true;
  for (auto i: range(D3Q15::NV)) {
    ans &= D3Q15::CX[i] == D3Q15::C[i][0];
    ans &= D3Q15::CY[i] == D3Q15::C[i][1];
    ans &= D3Q15::CZ[i] == D3Q15::C[i][2];
  }
  return ans;
}
static_assert(check_vels());

// Distribution function
static_assert(detail::is_component<DistFn>::value);
static_assert(field_traits<DistFn>::NELEM == 15);

// Density derived var
static_assert(detail::blah(static_cast<Density*>(nullptr)));
static_assert(detail::is_derived_var_v<Density>);
static_assert(field_traits<Density>::NELEM == 1);

int main() {
  // Just return success - if this executable gets compiled i
  return 0;
}
