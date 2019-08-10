// -*- mode: c++; -*-
#ifndef LAMBREX_VELOCITY_SET_H
#define LAMBREX_VELOCITY_SET_H

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

template <typename Impl, int ND_, int NV_, int HALO_>
struct VelocitySet {
  static constexpr int ND = ND_;
  static constexpr int NV = NV_;
  static constexpr int HALO = HALO_;
  using idx_type = std::array<int, ND>;
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
    std::array<idx_type, NV> ans = {0};
    for (auto i = 0; i < NV; ++i)
      ans[i] = {derived::CX[i], derived::CY[i],derived::CZ[i]};
    return ans;
  }
  static constexpr std::array<std::array<double, ND>, NV> make_vel_float_array() {
    std::array<std::array<double, ND>, NV> ans = {0};
    for (auto i = 0; i < NV; ++i)
      ans[i] = {derived::CX[i], derived::CY[i], derived::CZ[i]};
    return ans;
  }
  
public:
  static constexpr std::array<std::array<int, ND>, NV> CI = make_vel_int_array();
  static constexpr std::array<std::array<double, ND>, NV> XI = make_vel_float_array();

  static constexpr indexer index_begin();
  static constexpr indexer index_end();

};

#endif
