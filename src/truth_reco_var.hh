#ifndef TRUTH_RECO_VAR_HH
#define TRUTH_RECO_VAR_HH

template <typename T>
struct var {
  T det, truth;

  // operator % applies function f to both values
  template <typename F>
  inline auto operator%(F f) const noexcept(noexcept(f(std::declval<T>())))
  -> var<decltype(f(std::declval<T>()))> {
#ifndef VAR_ALWAYS_MC
    if (is_mc) return { f(det), f(truth) };
    else return { f(det), { } };
#else
    return { f(det), f(truth) };
#endif
  }

#define VAR_OP(OP) \
  template <typename U> \
  inline auto operator OP (const U& x) const \
  noexcept(noexcept(std::declval<T>() OP std::declval<U>())) \
  -> var<decltype(std::declval<T>() OP std::declval<U>())> { \
    return { (det OP x), (truth OP x) }; } \
  template <typename U> \
  inline auto operator OP (const var<U>& x) const \
  noexcept(noexcept(std::declval<T>() OP std::declval<U>())) \
  -> var<decltype(std::declval<T>() OP std::declval<U>())> { \
    return { (det OP x.det), (truth OP x.truth) }; }

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

  // bool conversion operator only for var<bool>
  template <typename U = bool>
  inline operator std::enable_if_t<
    std::is_same<T,bool>::value && std::is_same<U,bool>::value,
    bool> () noexcept { return det; }
};

// overload abs
template <typename T>
inline var<T> abs(const var<T>& x) noexcept {
  return { std::abs(x.det), std::abs(x.truth) };
}

template <typename T>
struct var<TTreeReaderValue<T>> {
  TTreeReaderValue<T> det;
#ifndef VAR_ALWAYS_MC
  std::experimental::optional<TTreeReaderValue<T>> truth;
#else
  TTreeReaderValue<T> truth;
#endif
  var(TTreeReader& tr, const std::string& name)
  : det(tr,("HGamEventInfoAuxDyn."+name).c_str())
#ifdef VAR_ALWAYS_MC
    , truth(tr,("HGamTruthEventInfoAuxDyn."+name).c_str()) { }
#else
  {
    if (is_mc) truth.emplace(tr,("HGamTruthEventInfoAuxDyn."+name).c_str());
  }
#endif
  using ret_type = std::conditional_t<
    std::is_floating_point<T>::value, double, T>;
  inline var<ret_type> operator*() {
#ifndef VAR_ALWAYS_MC
    if (truth) return { *det, **truth };
    else return { *det, { } };
#else
    return { *det, *truth };
#endif
  }

  // operator % applies function f to both values
  template <typename F>
  inline auto operator%(F f) noexcept(noexcept(f(std::declval<T>())))
  -> var<decltype(f(std::declval<T>()))> {
#ifndef VAR_ALWAYS_MC
    if (truth) return { f(*det), f(**truth) };
    else return { f(*det), { } };
#else
    return { f(*det), f(*truth) };
#endif
  }

#ifndef VAR_ALWAYS_MC
#define VAR_OP(OP) \
  template <typename U> \
  inline auto operator OP (const U& x) \
  noexcept(noexcept(std::declval<T>() OP std::declval<U>())) \
  -> var<decltype(std::declval<T>() OP std::declval<U>())> { \
    if (truth) return { ((*det) OP x), ((**truth) OP x) }; \
    else return { ((*det) OP x), { } }; \
  } \
  template <typename U> \
  inline auto operator OP (const var<U>& x) \
  noexcept(noexcept(std::declval<T>() OP std::declval<U>())) \
  -> var<decltype(std::declval<T>() OP std::declval<U>())> { \
    if (truth) return { ((*det) OP x), ((**truth) OP x.truth) }; \
    else return { ((*det) OP x), { } }; \
  }
#else
#define VAR_OP(OP) \
  template <typename U> \
  inline auto operator OP (const U& x) \
  noexcept(noexcept(std::declval<T>() OP std::declval<U>())) \
  -> var<decltype(std::declval<T>() OP std::declval<U>())> { \
    return { ((*det) OP x), ((*truth) OP x) }; } \
  template <typename U> \
  inline auto operator OP (const var<U>& x) \
  noexcept(noexcept(std::declval<T>() OP std::declval<U>())) \
  -> var<decltype(std::declval<T>() OP std::declval<U>())> { \
    return { ((*det) OP x.det), ((*truth) OP x.truth) }; }
#endif

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
};

// overload abs
template <typename T>
inline var<T> abs(var<TTreeReaderValue<T>>& x) noexcept {
#ifndef VAR_ALWAYS_MC
    if (x.truth) return { std::abs(*x.det), std::abs(**x.truth) };
    else return { std::abs(*x.det), { } };
#else
  return { std::abs(*x.det), std::abs(*x.truth) };
#endif
}

#endif
