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
  };
  
  template <size_t... I>
  Accessor<I...> accessor_deducer(std::index_sequence<I...>);
}

using Accessor = decltype(detail::accessor_deducer(std::make_index_sequence<AMREX_SPACEDIM>{}));


// Call the supplied function for every point in every box in the
// MultiFab. The function will be passed an Accessor instance which
// will index into a MultiFab for you.
template<typename PosFunc>
void for_point_in(const amrex::MultiFab& mf, PosFunc&& f) {
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const auto box = mfi.validbox();
    const auto& lo = box.smallEnd();
    const auto& hi = box.bigEnd();

    auto acc = Accessor{mfi};
    for     (int k = lo[2]; k <= hi[2]; ++k)
      for   (int j = lo[1]; j <= hi[1]; ++j)
	for (int i = lo[0]; i <= hi[0]; ++i) {
	  acc.pos = {i, j, k};
	  f(acc);
	}
  }
}

#endif
