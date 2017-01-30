#ifndef IVANP_UTILITY_HH
#define IVANP_UTILITY_HH

#include <string>

namespace ivanp {

template <typename T>
struct named_ptr {
  using type = T;
  type *p;
  std::string name;

  named_ptr(): p(nullptr), name() { }
  named_ptr(const named_ptr& n) = default;
  named_ptr(named_ptr&& n) = default;
  template <typename P, typename N>
  named_ptr(P&& ptr, N&& name)
  : p(std::forward<P>(ptr)), name(std::forward<N>(name)) { }

  inline type& operator*() const noexcept { return *p; }
  inline type* operator->() const noexcept { return p; }
};

#ifdef _GLIBCXX_VECTOR
template <typename T>
auto reserve(size_t n) {
  std::vector<T> _v;
  _v.reserve(n);
  return _v;
}
#endif

template <typename B>
constexpr unsigned check(B first) noexcept { return first ? 1u : 0u; }
template <typename B, typename... BB>
constexpr std::enable_if_t<sizeof...(BB),unsigned>
check(B first, BB... checks) noexcept {
  return first ? 1u+check(checks...) : 0u;
}

#ifdef _GLIBCXX_TUPLE
template <typename... T, size_t... I>
auto subtuple(const std::tuple<T...>& t, std::index_sequence<I...>) {
  static_assert(sizeof...(I) <= sizeof...(T),"");
  return std::make_tuple( std::get<I>(t)... );
}

template <typename F, typename... T, size_t... I>
inline decltype(auto)
call(F f, const std::tuple<T...>& t, std::index_sequence<I...>) {
  static_assert(sizeof...(I) <= sizeof...(T),"");
  return f( std::get<I>(t)... );
}
#endif

}

#endif
