// -*- mode: c++; -*-
#ifndef LAMBREX_PACK_TOOLS_H
#define LAMBREX_PACK_TOOLS_H

// Helpers for dealing with parameter packs


// Turn a pack into a single type
template <typename... Pack>
struct pack_wrapper {
};



// Map from type (assuming all T in Pack are unique) to its index
// within a pack. Result will be in a static member `value'.

// Primary template
template <typename T, typename... Pack>
struct t2i;

// Specialisation for when the target is at the front and have zero or
// more further types.
template <typename T, typename... Rest>
struct t2i<T, T, Rest...> {
  static constexpr std::size_t value = 0;
};

// Specialisation for target not at front
template <typename T, typename P0, typename... Rest>
struct t2i<T, P0, Rest...> {
  static constexpr std::size_t value = 1 + t2i<T, Rest...>::value;
};

// Error case
template <typename T>
struct t2i<T> {
  static_assert(sizeof(T) == 0, "Type not in pack");
};

// Case when we pass a pack wrapper
template<typename T, typename... Pack>
struct t2i<T, pack_wrapper<Pack...>> {
  static constexpr std::size_t value = t2i<T, Pack...>::value;
};

// Convenience variable template
template <typename T, typename... Pack>
inline constexpr std::size_t t2i_v = t2i<T, Pack...>::value;


// Get type by index - result in member typedef `type'

// Primary template
template <size_t i, typename... Pack>
struct i2t;

// Terminating case
template<typename Head, typename... Rest>
struct i2t<0, Head, Rest...> {
  using type = Head;
};
// Recursive case
template<size_t i, typename Head, typename... Rest>
struct i2t<i, Head, Rest...> {
  using type = typename i2t<i-1, Rest...>::type;
};
// pack wrapper
template<size_t i, typename... Pack>
struct i2t<i, pack_wrapper<Pack...>> {
  using type = typename i2t<i, Pack...>::type;
};

// Convenience alias
template <size_t i, typename... Pack>
using i2t_t = typename i2t<i, Pack...>::type;

#endif
