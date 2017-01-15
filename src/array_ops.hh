#ifndef IVANP_ARRAY_OPS
#define IVANP_ARRAY_OPS

#include <array>

template <typename T, size_t N, typename X>
inline std::array<T,N>& operator*=(std::array<T,N>& a, X x)
noexcept(noexcept(std::get<0>(a) *= std::declval<X>()))
{
  for (size_t i=0; i<N; ++i) a[i] *= x;
  return a;
}

template <typename T, typename X>
inline bool in(X x, const std::array<T,2>& a)
noexcept(
  noexcept(std::declval<T>() < std::declval<X>()) &&
  noexcept(std::declval<X>() < std::declval<T>()) )
{
  return (a[0] < x && x < a[1]);
}

template <typename T>
inline T len(const std::array<T,2>& a) noexcept {
  return (std::get<1>(a) - std::get<0>(a));
}

#endif
