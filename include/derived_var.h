// -*- mode: c++; -*-
#ifndef LAMBREX_DERIVED_VAR_H
#define LAMBREX_DERIVED_VAR_H

#include <AMReX_Box.H>
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

#include "pack_tools.h"
#include "field.h"

struct ScalarTag {
  static constexpr std::size_t RANK = 0;
  static constexpr std::size_t SIZE = 1;
};
  
// template <std::size_t DIMS...>
// struct ElementTag {
//   using pack_t = pack_wrapper<DIMS...>;
  
//   static constexpr std::size_t RANK = sizeof...(DIMS);
//   static constexpr std::size_t SIZE = (1 * ... * DIMS);
//   static_assert(SIZE > 0, "Invalid dimension of zero");
// };

// using ScalarTag = ElementTag<>;

// template <std::size_t DIM>
// using VectorTag = ElementTag<D>;

// CRTP base for "Derived Variables"
//
// A DV is something that is computed from one or more fields, whose
// types are passed to this as the dependencies.
//
// The implementation is the first parameter and must provide a member
// function (template)
// 
// void Impl::calculate(int x, int y, int z,
//                      SOMETHING... sources,
//                      elem_ref ans)
// 
// that will compute the DV at the given (x,y,z) point from the fields
// it depends on (sources...) in order of the tag types in
// dependencies. The result should be stored in ans(x,y,z, element_inds...)
// 
// This is followed by a tag class for the dimensionality of the
// computed quantity probably an intantiation of ElementTag
// (ScalarTag, VectorTag<N>, etc)
//
// Dependencies are zero or more fields that are inputs.
//
// TODO: specify the halo size for each dependency. Can't have 2
// parameter packs so prob have to wrap each dep in a tag template?
template <typename Impl, typename DIM, typename... dependencies>
struct DerivedVar {
  // static_assert(
  //   std::is_invocable_r<
  //     void,
  //     decltype(&Impl::calculate), Impl,
  //   dependencies::const_array_ref...>::value,
  //   "Implementation must have suitable calculate function");
  using dim_type = DIM;
  using pack_type = pack_wrapper<dependencies...>;
  using base_type = DerivedVar<Impl, DIM, dependencies...>;
  static constexpr std::size_t HALO = 0;

private:
  // CRTP helpers
  using derived = Impl;
  derived& self() {
    return *static_cast<derived*>(this);
  }
  const derived& self() const {
    return *static_cast<const derived*>(this);
  }
public:
  // Given a box defining the space, the list of source data and the
  // output array, do the computations.
  void fill_box(const amrex::Box& box,
		const typename field_traits<dependencies>::box_type&... sources,
		amrex::Array4<amrex::Real>& dest) {
  //void fill_box(const amrex::Box& box, typename dependencies::const_array_ref... sources, amrex::Array4<amrex::Real>& dest) {
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for     (int k = lo.z; k <= hi.z; ++k)
      for   (int j = lo.y; j <= hi.y; ++j)
	for (int i = lo.x; i <= hi.x; ++i)
	  self().calculate(i, j, k, sources..., dest);
  }
};


// template <typename T, int NC>
// struct ElemWrapper {
//   amrex::Array4<T> data;
//   int i, j, k;

//   std::enable_if_t< (NC>1), T&> operator[](int c) {
//     return data(i,j,k, c);
//   }

//   std::enable_if_t< NC==1, ElemWrapper&> operator=(T& x) {
//     data(i,j,k, 0) = x;
//     return *this;
//   }
//   std::enable_if_t< NC==1, ElemWrapper&> operator+=(T& x) {
//     data(i,j,k, 0) += x;
//     return *this;
//   }
// };

#endif
