// -*- mode: c++; -*-
#ifndef LAMBREX_COMPONENT_H
#define LAMBREX_COMPONENT_H

template <typename VS>
struct Component {
  using VelocitySet = VS;
  static constexpr auto ND = VS::ND;
  static constexpr auto NV = VS::NV;
  static constexpr auto CI = VS::CI;
  static constexpr auto XI = VS::XI;
  static constexpr auto HALO = VS::HALO;

  static void DefineLevelData(amrex::MultiFab& fab, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
    fab.define(ba, dm, NV, HALO);
  }
  static amrex::MultiFab MakeLevelData(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
    return amrex::MultiFab{ba, dm, NV, HALO};
  }
};


#endif
