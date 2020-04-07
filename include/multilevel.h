// -*- mode: c++; -*-
#ifndef LAMBREX_MULTILEVEL_H
#define LAMBREX_MULTILEVEL_H

#include <utility>
#include <array>
#include <algorithm>

#include <AMReX_REAL.H>
#include <AMReX_MultiFab.H>

#include "pack_tools.h"
#include "field.h"

// A single level's time info
struct TimeData {
  amrex::Real current;
  amrex::Real delta;
  int step;
};

// Store a single level's fields at a single time instant in a
// sequence of MultiFabs
template <typename... Fields>
struct State {
  static constexpr auto NFIELDS = sizeof...(Fields);
  using index = std::index_sequence_for<Fields...>;
  using field_types = pack_wrapper<Fields...>;

  using MFArray = std::array<amrex::MultiFab, NFIELDS>;
  MFArray fields;

  template <typename F>
  amrex::MultiFab& get() {
    constexpr auto i = t2i_v<F, Fields...>;
    return fields[i];
  }
  template <typename F>
  const amrex::MultiFab& get() const {
    constexpr auto i = t2i_v<F, Fields...>;
    return fields[i];
  }

  void Define(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
    auto x = {(Def<Fields>(ba, dm), 0)...};
  }

  template<typename F>
  void Def(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
    constexpr auto i = t2i_v<F, Fields...>;
    field_traits<F>::DefineLevelData(fields[i], ba, dm);
  }
  void Clear() {
    for (auto& f: fields) {
      f.clear();
    }
  }
};

// Store the state at a given level at the current time instant and
// provide a target for the next time instant's state too
template <typename State>
struct LevelData {
  void Define(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
    now.Define(ba, dm);
    next.Define(ba, dm);
  }
  void Clear() {
    now.Clear();
    next.Clear();
    time.current = 0.0;
    time.delta = 0.0;
    time.step = 0;
  }
  void UpdateNow() {
    // now becomes next
    // will need to be more sophisticated when there are more than two steps
    std::swap(now, next);
    return;
  }

  TimeData time;
  State now;
  State next;
};


#endif
