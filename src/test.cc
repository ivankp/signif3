#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <tuple>
#include <memory>
#include <regex>

#include "binner.hh"
#include "re_axes.hh"
#include "catstr.hh"
#include "timed_counter.hh"
#include "array_ops.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#include "migration_axis.hh"

using std::cout;
using std::cerr;
using std::endl;

template <typename T>
struct var {
  T det, truth;

#define VAR_OP(OP) \
  inline auto operator OP (const T& x) const noexcept(noexcept(det OP x)) \
  -> var<decltype(det OP x)> { return { (det OP x), (truth OP x) }; }

  VAR_OP(==)
  VAR_OP(!=)
  VAR_OP(<)
  VAR_OP(<=)
  VAR_OP(>)
  VAR_OP(>=)

  template <typename F>
  inline auto operator()(F f) const noexcept(noexcept(f(std::declval<T>())))
  -> var<decltype(f(std::declval<T>()))> { return { f(det), f(truth) }; }
};

using hist = ivanp::binner<double, ivanp::tuple_of_same_t<
  // ivanp::axis_spec<migration_axis<double>,false,false>,2>>;
  ivanp::axis_spec<migration_axis<double>,true,true>,2>>;

auto make_hist(const char* name, const re_axes& ra, unsigned nchecks) {
  migration_axis<double> axis(&*ra[name],nchecks);
  return hist( name, axis, axis );
}

template <typename T, typename... B>
void fill(hist& h, double w, const var<T>& x, const var<B>&... checks) {
  auto hbin = 
  h( std::tie(x.det,   checks.det...),
     std::tie(x.truth, checks.truth...), w );
  test( hbin )
  std::cout << std::endl;
}

int main(int argc, char* argv[])
{
  cout << std::boolalpha;
  re_axes ra("test.bins");

  auto h = make_hist("test",ra,1);

  fill(h,1,var<double>{1,1},var<bool>{true,true});
  fill(h,2,var<double>{2,2},var<bool>{false,true});
  fill(h,5,var<double>{3,3},var<bool>{true,false});
  fill(h,9,var<double>{4,4},var<bool>{false,false});
  fill(h,13,var<double>{5,20},var<bool>{true,true});
  fill(h,17,var<double>{20,6},var<bool>{true,true});
  fill(h,50,var<double>{20,20},var<bool>{true,true});

  for (const auto& h : hist::all) {
    cout << h.name << endl;
    for (auto& b : h->bins()) cout << b << endl;
    cout << endl;
  }

  return 0;
}


