#ifndef TRUTH_RECO_VAR_HH
#define TRUTH_RECO_VAR_HH

#ifndef VAR_ALWAYS_MC
#include <experimental/optional>
#endif

#include "tuple_comprehension.hh"
#include "catstr.hh"

template <typename T>
struct var {
  T det, truth;

  // operator | applies function f to both values
  template <typename F>
  inline auto operator|(F&& f) const noexcept(noexcept(f(det)))
  -> var<decltype(f(det))> {
    return { f(det),
#ifndef VAR_ALWAYS_MC
      is_mc ? f(truth) : decltype(f(truth)){}
#else
      f(truth)
#endif
    };
  }

  template <typename F1, typename F2>
  inline auto operator|(const std::pair<F1,F2>& f) const
  noexcept(noexcept(f.first(det)) && noexcept(f.second(truth)))
  -> var<decltype(f.first(det))> {
    return { f.first(det),
#ifndef VAR_ALWAYS_MC
      is_mc ? f.second(truth) : decltype(f.second(truth)){}
#else
      f.second(truth)
#endif
    };
  }

#define VAR_OP(OP) \
  template <typename U> \
  inline auto operator OP (const U& x) const \
  noexcept(noexcept(det OP x)) -> var<decltype(det OP x)> { \
    return { det OP x, truth OP x }; \
  } \
  template <typename U> \
  inline auto operator OP (const var<U>& x) const \
  noexcept(noexcept(det OP x.det)) -> var<decltype(det OP x.det)> { \
    return { det OP x.det, truth OP x.truth }; \
  }

  VAR_OP(==)
  VAR_OP(!=)
  VAR_OP(<)
  VAR_OP(<=)
  VAR_OP(>)
  VAR_OP(>=)
  VAR_OP(+)
  VAR_OP(-)
  VAR_OP(*)
  VAR_OP(/)

#undef VAR_OP

#define VAR_OP(OP) \
  template <typename U> \
  inline var& operator OP (const U& x) \
  noexcept(noexcept(det OP x)) { \
    det OP x; truth OP x; return *this; \
  } \
  template <typename U> \
  inline var& operator OP (const var<U>& x) \
  noexcept(noexcept(det OP x.det)) { \
    det OP x.det; truth OP x.truth; return *this; \
  }

  VAR_OP(+=)
  VAR_OP(-=)
  VAR_OP(*=)
  VAR_OP(/=)

#undef VAR_OP

  // bool conversion operator only for var<bool>
  template <typename U = std::decay_t<T>>
  inline operator std::enable_if_t< std::is_same<U,bool>::value,
  bool> () noexcept { return det; }

  // overload abs
  friend inline var abs(const var& x) noexcept {
    return { std::abs(x.det), std::abs(x.truth) };
  }
};

#ifdef ROOT_TTreeReaderValue

template <typename T>
class var<TTreeReaderValue<T>> {
  TTreeReaderValue<T> _det;
#ifndef VAR_ALWAYS_MC
  std::experimental::optional<TTreeReaderValue<T>> _truth;
#else
  TTreeReaderValue<T> _truth;
#endif

public:
  var(TTreeReader& tr, const std::string& name)
  : _det(tr,("HGamEventInfoAuxDyn."+name).c_str())
#ifdef VAR_ALWAYS_MC
    , _truth(tr,("HGamTruthEventInfoAuxDyn."+name).c_str()) { }
#else
  {
    if (is_mc) _truth.emplace(tr,("HGamTruthEventInfoAuxDyn."+name).c_str());
  }
#endif

  using type = std::conditional_t<
    std::is_floating_point<T>::value, double, T>;

  inline type det() { return *_det; }
#ifndef VAR_ALWAYS_MC
  inline type truth() {
    if (_truth) return **_truth;
    else return { };
  }
#else
  inline type truth() { return *_truth; }
#endif

  inline var<type> operator*() { return { det(), truth() }; }

  // operator | applies function f to both values
  template <typename F>
  inline auto operator|(F&& f) noexcept(noexcept(f(std::declval<type>())))
  -> var<decltype(f(std::declval<type>()))> { return { det(), truth() }; }

#define VAR_OP(OP) \
  template <typename U> \
  inline auto operator OP (const U& x) \
  noexcept(noexcept(std::declval<type>() OP x)) \
  -> var<decltype(std::declval<type>() OP x)> { \
    return { det() OP x, truth() OP x }; \
  }

  VAR_OP(==)
  VAR_OP(!=)
  VAR_OP(<)
  VAR_OP(<=)
  VAR_OP(>)
  VAR_OP(>=)
  VAR_OP(+)
  VAR_OP(-)
  VAR_OP(*)
  VAR_OP(/)

#undef VAR_OP

  // overload abs
  friend inline var<type> abs(var& x) noexcept {
    return { std::abs(x.det()), std::abs(x.truth()) };
  }
};

#endif

#ifdef ROOT_TTreeReaderArray

template <typename T, size_t N>
class var<std::array<TTreeReaderArray<T>,N>> {
  using reader = TTreeReaderArray<T>;
  template <typename U> using array = std::array<U,N>;

  array<reader> _det;
#ifndef VAR_ALWAYS_MC
  std::experimental::optional<array<reader>>
#else
  array<reader>
#endif
    _truth;

#define MAKE_READER(I) \
  [&](const std::string& name) -> reader { \
    return { tr, (std::get<I>(obj)+name).c_str() }; \
  }

public:
  var(TTreeReader& tr,
      const std::array<std::string,2>& obj,
      const array<std::string>& names)
  : _det(names | MAKE_READER(0))
#ifdef VAR_ALWAYS_MC
  , _truth(names | MAKE_READER(1)) { }
#else
  {
    if (is_mc) _truth.emplace( names | MAKE_READER(1) );
  }
#endif

  var(TTreeReader& tr,
      const std::array<std::string,2>& obj,
      const array<std::string>& names_det,
      const array<std::string>& names_truth)
  : _det(names_det | MAKE_READER(0))
#ifdef VAR_ALWAYS_MC
  , _truth(names_truth | MAKE_READER(1)) { }
#else
  {
    if (is_mc) _truth.emplace( names_truth | MAKE_READER(1) );
  }
#endif

#undef MAKE_READER

  using type = std::conditional_t<
    std::is_floating_point<T>::value, double, T>;

  inline var<array<type>> operator[](unsigned i) {
    const auto f = [i](reader& x) -> type {
      // if (i >= x.GetSize()) throw std::runtime_error(ivanp::cat(
      //   x.GetBranchName()," index (",i,") out of range (",x.GetSize(),')'));
      if (i >= x.GetSize()) return { };
      return x[i];
    };
    return { _det | f,
#ifndef VAR_ALWAYS_MC
      _truth ? (*_truth) | f : array<type>{ }
#else
      _truth | f
#endif
    };
  }
};

#endif

#endif
