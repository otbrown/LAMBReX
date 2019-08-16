// -*- mode: c++; -*-
#ifndef LAMBREX_VELOCITY_SET_H
#define LAMBREX_VELOCITY_SET_H

#include <array>
#include <cmath>

#include "range.hpp"

template <typename Int>
class IdxIterator {
  Int i = 0;

public:
  constexpr IdxIterator() {
  }
  constexpr IdxIterator(Int start) : i(start) {
  }
  constexpr IdxIterator& operator++() {
    ++i;
    return *this;
  }
  constexpr Int operator*() const {
    return i;
  }
};

template <typename Int>
constexpr bool operator==(const IdxIterator<Int>& a, const IdxIterator<Int>& b) {
  return a.i == b.i;
}
template <typename Int>
constexpr bool operator!=(const IdxIterator<Int>& a, const IdxIterator<Int>& b) {
  return !(a == b);
}
template <typename Int>
constexpr bool operator<(const IdxIterator<Int>& a, const IdxIterator<Int>& b) {
  return a.i < b.i;
}
template <typename Int>
constexpr bool operator>(const IdxIterator<Int>& a, const IdxIterator<Int>& b) {
  return b < a;
}
template <typename Int>
constexpr bool operator<=(const IdxIterator<Int>& a, const IdxIterator<Int>& b) {
  return !(a > b);
}
template <typename Int>
constexpr bool operator>=(const IdxIterator<Int>& a, const IdxIterator<Int>& b) {
  return !(a < b);
}

constexpr int kronecker(int a, int b) {
  return a == b ? 1 : 0;
}

template <typename Impl, int ND_, int NV_, int HALO_>
struct VelocitySet {
  static constexpr int ND = ND_;
  static constexpr int NV = NV_;
  static constexpr int HALO = HALO_;
  static constexpr double CS2 = 1.0 / 3.0;
  static constexpr double CS = std::sqrt(CS2);

  using idx_type = std::array<int, ND>;
  using vector_type = std::array<double, ND>;
  using matrix_type = std::array<vector_type, ND>;

  using indexer = IdxIterator<int>;
private:
  using derived = Impl;

  derived& self() {
    return *static_cast<derived*>(this);
  }
  const derived& self() const {
    return *static_cast<const derived*>(this);
  }

  static constexpr std::array<idx_type, NV> make_vel_int_array() {
    std::array<idx_type, NV> ans{};
    for (auto i: range(NV))
      ans[i] = {derived::CX[i], derived::CY[i],derived::CZ[i]};
    return ans;
  }
  static constexpr auto make_vel_float_array() {
    std::array<vector_type, NV> ans{};
    for (auto i: range(NV))
      ans[i] = {derived::CX[i], derived::CY[i], derived::CZ[i]};
    return ans;
  }
  static constexpr auto make_q_array() {
    std::array<matrix_type, NV> ans{};
    for (auto i: range(NV))
      for (auto a: range(ND))
	for (auto b: range(ND))
	  ans[i][a][b] = X[i][a]*X[i][b] - CS2 * kronecker(a,b);
    return ans;
  }
public:
  static constexpr std::array<idx_type, NV> C = make_vel_int_array();
  static constexpr std::array<vector_type, NV> X = make_vel_float_array();
  static constexpr std::array<matrix_type, NV> Q = make_q_array();

  static constexpr indexer index_begin();
  static constexpr indexer index_end();

};

#endif
