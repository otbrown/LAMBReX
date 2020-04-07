// -*- mode: c++; -*-
#ifndef LAMBREX_AMR_HELP_H
#define LAMBREX_AMR_HELP_H

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_IntVect.H>

namespace detail {
  // This type will index into a MultiFab at the desired point and
  // also provides shifts (by operators + and -) for stencil
  // operations.
  //
  // For this to work, the pack I must be an index_sequence up to
  // AMREX_SPACEDIM, i.e. DIM=2 => I={0,1} and DIM=3 => I={0,1,2}
  //
  // We use the undefined deducer function below to produce the
  // correct alias for us in the main namespace.
  //
  // Example:
  // MultiFab f, rho;
  // for (auto mf_iter = ...){
  //   for (x,y,z = ...) {
  //     auto point = Accessor{mf_iter, {x,y,z}};
  //     point(f, i) = w[i]* point(rho);
  //   }
  // }
  template<size_t... I>
  struct Accessor {
    amrex::MFIter& mfi;
    amrex::IntVect pos;

    inline double& operator()(amrex::MultiFab& multifab, int comp=0) const {
      return multifab[mfi](pos, comp);
    }
    inline const double& operator()(const amrex::MultiFab& multifab, int comp=0) const {
      return multifab[mfi](pos, comp);
    }

    inline int& operator()(amrex::iMultiFab& imultifab, int comp=0) const {
      return imultifab[mfi](pos, comp);
    }
    inline const int& operator()(const amrex::iMultiFab& imultifab, int comp=0) const {
      return imultifab[mfi](pos, comp);
    }

    Accessor operator+(const amrex::IntVect& shift) const {
      return Accessor{mfi, pos+shift};
    }
    Accessor operator-(const amrex::IntVect& shift) const {
      return Accessor{mfi, pos-shift};
    }
    inline Accessor operator+(const std::array<int, AMREX_SPACEDIM>& shift) const {
      return Accessor{mfi, {(pos[I] + shift[I])...}};
    }
    Accessor operator-(const std::array<int, AMREX_SPACEDIM>& shift) const {
      return Accessor{mfi, {(pos[I] - shift[I])...}};
    }

    Accessor CoarseAccess(const int RATIO) const {
      amrex::IntVect coarse_pos(pos);
      coarse_pos.coarsen(RATIO);

      return Accessor{mfi, coarse_pos};
    }
  };

  template <size_t... I>
  Accessor<I...> accessor_deducer(std::index_sequence<I...>);
}

using Accessor = decltype(detail::accessor_deducer(std::make_index_sequence<AMREX_SPACEDIM>{}));


// Call the supplied function for every point in every box in the
// MultiFab. The function will be passed an Accessor instance which
// will index into a MultiFab for you.
template<typename PosFunc>
void for_point_in(const amrex::MultiFab& mf, PosFunc&& f, int const HALO_DEPTH = 0) {
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const auto box = mfi.validbox();
    const auto& lo = box.smallEnd();
    const auto& hi = box.bigEnd();

    auto acc = Accessor{mfi};
    // loops ordered fortran style
    for (int i = lo[0]-HALO_DEPTH; i <= hi[0]+HALO_DEPTH; ++i)
      for (int j = lo[1]-HALO_DEPTH; j <= hi[1]+HALO_DEPTH; ++j)
        for (int k = lo[2]-HALO_DEPTH; k <= hi[2]+HALO_DEPTH; ++k) {
	         acc.pos = {i, j, k};
	         f(acc);
         }
  }
}

// Call the supplied function for every point in the boundary of every box in
// the MultiFab. The function will be passed an Accessor instance which
// will index into a MultiFab for you.
template<typename PosFunc>
void for_point_in_boundary(const amrex::MultiFab& mf, int const HALO_DEPTH,
PosFunc&& f) {
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    int i, j, k;
    const auto box = mf[mfi].box();
    const auto& lo = box.smallEnd();
    const auto& hi = box.bigEnd();

    auto acc = Accessor{mfi};
    // loops ordered fortran style
    // lower j-k planes
    for (i = lo[0]; i < lo[0]+HALO_DEPTH; ++i) {
      for (j = lo[1]; j <= hi[1]; ++j) {
        for (k = lo[2]; k <= hi[2]; ++k) {
          acc.pos = {i, j, k};
          f(acc);
        }
      }
    }
    // upper j-k planes
    for (i = hi[0]-HALO_DEPTH+1; i <= hi[0]; ++i) {
      for (j = lo[1]; j <= hi[1]; ++j) {
        for (k = lo[2]; k <= hi[2]; ++k) {
          acc.pos = {i, j, k};
          f(acc);
        }
      }
    }

    // lower i-k planes
    for (j = lo[1]; j < lo[1]+HALO_DEPTH; ++j) {
      for (i = lo[0]+HALO_DEPTH; i <= hi[0]-HALO_DEPTH; ++i) {
        for (k = lo[2]; k <= hi[2]; ++k) {
          acc.pos = {i, j, k};
          f(acc);
        }
      }
    }
    // upper i-k planes
    for (j = hi[1]-HALO_DEPTH+1; j <= hi[1]; ++j) {
      for (i = lo[0]+HALO_DEPTH; i <= hi[0]-HALO_DEPTH; ++i) {
        for (k = lo[2]; k <= hi[2]; ++k) {
          acc.pos = {i, j, k};
          f(acc);
        }
      }
    }

    // lower i-j planes
    for (k = lo[2]; k < lo[2]+HALO_DEPTH; ++k) {
      for (i = lo[0]+HALO_DEPTH; i <= hi[0]-HALO_DEPTH; ++i) {
        for (j = lo[1]+HALO_DEPTH; j <= hi[1]-HALO_DEPTH; ++j) {
          acc.pos = {i, j, k};
          f(acc);
        }
      }
    }
    // upper i-j planes
    for (k = hi[2]-HALO_DEPTH+1; k <= hi[2]; ++k) {
      for (i = lo[0]+HALO_DEPTH; i <= hi[0]-HALO_DEPTH; ++i) {
        for (j = lo[1]+HALO_DEPTH; j <= hi[1]-HALO_DEPTH; ++j) {
          acc.pos = {i, j, k};
          f(acc);
        }
      }
    }
  }
}

#endif
