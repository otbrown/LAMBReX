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
template <typename VS>
struct Component {
  using VelocitySet = VS;
  static constexpr auto ND = VS::ND;
  static constexpr auto NV = VS::NV;
  static constexpr auto HALO = VS::HALO;
};

#endif
