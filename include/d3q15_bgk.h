// -*- mode: c++; -*-
#ifndef LAMBREX_D3Q15_BGK_H
#define LAMBREX_D3Q15_BGK_H

#include "velocity_set.h"
#include "component.h"
#include "derived_var.h"
#include "multilevel.h"
#include "range.hpp"

struct D3Q15 : public VelocitySet<D3Q15, 3, 15, 2> {
  using base = VelocitySet<D3Q15, 3, 15, 2>;

  static constexpr std::array<int, NV> CX = {
    0,    1, -1, 0, 0, 0, 0,    1, 1, 1, 1, -1, -1, -1, -1
  };
  static constexpr std::array<int, NV> CY = {
    0,    0, 0, 1, -1, 0, 0,    1, 1, -1, -1, 1, 1, -1, -1
  };
  static constexpr std::array<int, NV> CZ = {
    0,    0, 0, 0, 0, 1, -1,    1, -1, 1, -1, 1, -1, 1, -1
  };

  static constexpr std::array<double, NV> W = {
    2./9.,
    1./9., 1./9., 1./9., 1./9., 1./9., 1./9.,
    1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.
  };
};

struct DistFn : public Component<D3Q15> {
};

struct Density : DerivedVar<Density, ScalarTag, DistFn> {
  template <typename F, typename RHO>
  static void calculate(int x, int y, int z, F f, RHO ans) {
    ans = 0.0;
    for (auto i: range(DistFn::NV))
      ans += f(x, y, z, i);
  }
};

// template <typename DistT>
// struct momentum_density : DerivedVar<momentum_density, Vector<DistT::ND>, DistT> {
//   template <typename elem_ref>
//   void calculate(int x, int y, int z, const typename DistT::ArrayT& f, elem_ref ans) {
//     using VS = typename DistT::VS;
//     for (int alpha = 0; alpha < VS::ND; ++alpha) {
//       ans[alpha] = 0.0;
//       for (int i = 0; i < DistT::NV; ++i) {
// 	ans[alpha] += f(x, y, z, i) * VS::XI[i][alpha];
//       }
//     }
//   }
// };

using SimState = State<DistFn, Density>;
using SimLevelData = LevelData<SimState>;


#endif
