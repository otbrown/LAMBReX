// -*- mode: c++; -*-
#ifndef LAMBREX_FIELD_H
#define LAMBREX_FIELD_H

#include <type_traits>

// Forward declarations
#include <AMReX_MultiFab.H>
namespace amrex {
  //class MultiFab;
  class BoxArray;
  class DistributionMapping;
}

template <typename VS>
struct Component;
template <typename Impl, typename DIM, typename... dependencies>
struct DerivedVar;

// A field is a quanity that exists across space at a moment in
// time. It may store its values or not, but if it does, this traits
// class will explain how to do so to AMReX.

// Member types:
//
// field_type - the type passed in
//
// field_category - a tag type to identify what sort of field it is



// Static member variables:
//
// size_t NELEM - how many Reals are stored at each point in space
// 
// size_t HALO - the size of halo needed

// static member functions:
// 
// void DefineLevelData(amrex::MultiFab& leveldata, const amrex::BoxArray&, const amrex::DistributionMapping&);
// Initialise a given refinement level's MultiFab
//
// amrex::MultiFab MakeLevelData(const amrex::BoxArray&, const amrex::DistributionMapping&);
// Create a new multifab

template <typename Field>
struct field_traits;

namespace detail {
  // Tag types for category of field
  struct component_tag {
  };
  struct derived_var_tag {
  };

  // Trait for checking if a type is a Component
  template <typename T, typename Enable = void>
  struct is_component : public std::false_type {};
  template <typename T>
  struct is_component<T,
		      std::enable_if_t<std::is_base_of_v<Component<typename T::VelocitySet>, T>>
		      > : public std::true_type {};
  template <typename T>
  constexpr bool is_component_v = is_component<T>::value;

  template <typename Impl, typename DIM, typename... dependencies>
  constexpr bool blah(const DerivedVar<Impl, DIM, dependencies...>* ) {
    return true;
  }

  constexpr bool foo(...) {
    return false;
  }
  template <typename T>
  constexpr auto foo(const T* t = nullptr) -> decltype(blah(t)) {
    return true;
  }

  // // Trait for checking if a type is a DerivedVar
  // template <typename T, typename Enable = void>
  // struct is_derived_var : public std::false_type {};
  // template <typename T>
  // struct is_derived_var<T,
  // 			std::enable_if_t<std::is_base_of_v<DerivedVar<T, typename T::dim_type,  typename T::pack_type>, T>>
  // 		      > : public std::true_type {};

  template <typename T>
  constexpr bool is_derived_var_v = foo(static_cast<const T*>(nullptr));


  // Traits that have to be specialised live here
  //
  // Primary template
  template <typename Field, typename Enable=void>
  struct field_traits_impl;
  
  // Specialisation for Components
  template <typename Field>
  struct field_traits_impl<Field, std::enable_if_t<is_component_v<Field>>> {
    using field_type = Field;
    using field_category = component_tag;
    using box_type = amrex::Array4<amrex::Real>;
    
    static constexpr auto NELEM = Field::NV;
    static constexpr auto HALO = Field::HALO;
  };

  // Specialisation for DerivedVars
  template <typename Field>
  struct field_traits_impl<Field, std::enable_if_t<is_derived_var_v<Field>>> {
    using field_type = Field;
    using field_category = component_tag;
    using box_type = amrex::Array4<amrex::Real>;
    
    static constexpr auto NELEM = Field::dim_type::SIZE;
    static constexpr auto HALO = Field::HALO;
  };
}

template <typename Field>
struct field_traits : public detail::field_traits_impl<Field> {
  using base_type = detail::field_traits_impl<Field>;
  
  static void DefineLevelData(amrex::MultiFab& fab, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
    fab.define(ba, dm, field_traits::NELEM, field_traits::HALO);
  }
  static amrex::MultiFab MakeLevelData(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
    return amrex::MultiFab{ba, dm, field_traits::NELEM, field_traits::HALO};
  }
};

#endif
