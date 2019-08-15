// -*- mode: c++; -*-
#ifndef LAMBREX_COMPONENT_H
#define LAMBREX_COMPONENT_H

// Base for tag classes for a lattice Boltzmann component
//
// Create a subclass of an instantiation of this for every componen in
// your system.
//
// This class just aliases the VelocitySet type and uses it's `ND`,
// `NV`, and `HALO` members.
//
// This class template also has a generic implementation of
// propagation/streaming.
template <typename VS>
struct Component {
  using VelocitySet = VS;
  static constexpr auto ND = VS::ND;
  static constexpr auto NV = VS::NV;
  static constexpr auto HALO = VS::HALO;

  template<typename Accessor, typename MFin, typename MFout>
  static void PropagatePoint(const Accessor& dest_pt, const MFin& f_pc, MFout& f_new) {
    for (auto i: range(NV)) {
      // Minus because pull-style propagation
      auto src_pt = dest_pt - VelocitySet::C[i];
      dest_pt(f_new, i) = src_pt(f_pc, i);
    }
  }
};

#endif
