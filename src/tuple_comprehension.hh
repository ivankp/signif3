#ifndef IVANP_TUPLE_COMPREHENSION
#define IVANP_TUPLE_COMPREHENSION

#include <utility>
#include <tuple>
#include <array>
#include "meta.hh"

namespace detail { namespace comprehension {

template <typename F, typename Tup, size_t... I>
inline auto tuple_comprehension(Tup&& t, F&& f, std::index_sequence<I...>)
-> std::tuple<decltype(f(std::get<I>(std::forward<Tup>(t))))...>
{
  return { f(std::get<I>(std::forward<Tup>(t)))... };
}

template <typename F, typename Arr, size_t... I>
inline auto array_comprehension(Arr&& a, F&& f, std::index_sequence<I...>)
-> std::array<decltype(f(std::get<0>(std::forward<Arr>(a)))), sizeof...(I)>
{
  return { f(std::get<I>(std::forward<Arr>(a)))... };
}

template <typename F, typename T>
using detect_pass_as_const = decltype(
  std::declval<F>()(std::declval<const T&>()) );

template <typename F, typename T>
using enable_nonconst_t = std::enable_if_t<
  !::ivanp::is_detected_v< detect_pass_as_const, F, T >
>;

}} // end namespace detail

template <typename F, typename... Ts>
inline auto operator|(const std::tuple<Ts...>& t, F&& f) {
  return ::detail::comprehension::tuple_comprehension(
    t, std::forward<F>(f), std::index_sequence_for<Ts...>{} );
}

template <typename F, typename... Ts>
inline auto operator|(std::tuple<Ts...>&& t, F&& f) {
  return ::detail::comprehension::tuple_comprehension(
    std::move(t), std::forward<F>(f), std::index_sequence_for<Ts...>{} );
}

template <typename F, typename T, size_t N,
          typename = ::detail::comprehension::enable_nonconst_t<F,T> >
inline auto operator|(std::array<T,N>& a, F&& f) {
  return ::detail::comprehension::array_comprehension(
    a, std::forward<F>(f), std::make_index_sequence<N>{} );
}

template <typename F, typename T, size_t N>
inline auto operator|(const std::array<T,N>& a, F&& f) {
  return ::detail::comprehension::array_comprehension(
    a, std::forward<F>(f), std::make_index_sequence<N>{} );
}

template <typename F, typename T, size_t N>
inline auto operator|(std::array<T,N>&& a, F&& f) {
  return ::detail::comprehension::array_comprehension(
    std::move(a), std::forward<F>(f), std::make_index_sequence<N>{} );
}

#endif
