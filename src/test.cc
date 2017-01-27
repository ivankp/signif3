#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <memory>

#include "binner.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

const bool is_mc = true;
const bool is_fiducial = true;
const double factor = 1.;
const double fw = 1.;
const double lumi = 1.;
const double n_all_inv = 1.;

struct hist_bin {
  static double weight;
  double
    tmp, tmp_truth,
    sig, sig_truth,
    bkg, signif,
    sb_purity, // signal-background purity
    tr_purity; // truth-reco purity

  hist_bin():
    tmp(0), tmp_truth(0),
    sig(0), sig_truth(0),
    bkg(0), signif(0),
    sb_purity(0), tr_purity(0) { }
  inline hist_bin& operator++() noexcept {
    tmp += weight;
    return *this;
  };
  inline hist_bin& operator()(bool matches_truth) noexcept {
    tmp += weight;
    if (matches_truth) tmp_truth += weight;
    return *this;
  }
  inline hist_bin& operator+=(const hist_bin& b) noexcept {
    sig += b.sig;
    sig_truth += b.sig_truth;
    bkg += b.bkg;
    return *this;
  }

  void merge() noexcept {
    if (is_mc) {
      sig += tmp*n_all_inv;
      sig_truth += tmp_truth*n_all_inv;
      tmp_truth = 0;
    } else {
      bkg += tmp;
    }
    tmp = 0;
  }
  void compute() noexcept {
    bkg *= factor;
    sig *= lumi;
    sig_truth *= lumi;
    const double sb = sig + bkg;
    signif = sig/std::sqrt(sb);
    sb_purity = sig/sb;
    tr_purity = sig_truth/sig;
  }
};
double hist_bin::weight;

std::ostream& operator<<(std::ostream& o, const hist_bin& b) {
  const auto prec = o.precision();
  const std::ios::fmtflags f( o.flags() );
  o << round(b.sig) << ' '
    << round(b.bkg) << ' '
    << std::setprecision(2) << std::fixed
    << b.signif << ' '
    << 100*b.sb_purity << "% "
    << 100*b.tr_purity << '%'
    << std::setprecision(prec);
  o.flags( f );
  return o;
}

template <typename... Axes>
using hist = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;

template <typename T>
struct var {
  T det, truth;
  inline bool operator==(const T& x) const noexcept { return (det == x); }
  inline bool operator< (const T& x) const noexcept { return (det <  x); }
  inline bool operator> (const T& x) const noexcept { return (det >  x); }
};

template <typename T1, typename T2, typename A1, typename A2>
void fill(hist<A1,A2>& h, const var<T1>& x1, const var<T2>& x2, bool match=true) {
  const auto bin_det = h.find_bin(x1.det,x2.det);
  test( bin_det )
  if (is_mc) {
    const auto bin_truth = h.find_bin(x1.truth,x2.truth);
    h.fill_bin(bin_det, is_fiducial && (bin_det == bin_truth) && match);

    test( std::get<0>(h.axes()).find_bin(x1.det) )
    test( std::get<1>(h.axes()).find_bin(x2.det) )

  } else h.fill_bin(bin_det);
}

int main(int argc, char* argv[])
{
  hist_bin::weight = 1;

  
  using hist2 = hist<
    ivanp::container_axis<std::vector<double>>,
    ivanp::container_axis<std::vector<double>> >;
  hist2
    h_Dphi_Dy_jj("Dphi_Dy_jj",{0.,M_PI_2,M_PI},{0.,2.,8.8});

  // var<double> dphi{0.1,0.1}, dy{3.,3.};
  for (double x : {-1.,1.,3.,9.}) {
    var<double> dphi{0.1,0.1}, dy{x,x};

    fill(h_Dphi_Dy_jj, dphi, dy);

    for (auto& bin : h_Dphi_Dy_jj.bins()) {
      cout << bin.tmp << ' ';
    }
    cout << endl;
  }
  cout << endl;

  for (auto& bin : h_Dphi_Dy_jj.bins()) {
    bin.merge();
    bin.compute();
  }
  for (auto& bin : h_Dphi_Dy_jj.bins()) {
    cout << bin.sig << ' ';
  }
  cout << endl;
  cout << endl;

  for (int j=0; j<4; ++j) {
    for (int i=0; i<4; ++i) {
      cout << ' ' << h_Dphi_Dy_jj.index(i,j);
    }
    cout << endl;
  }
  cout << endl;

  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      cout << ' ' << h_Dphi_Dy_jj.index(i,j);
    }
    cout << endl;
  }
  cout << endl;

  for (int i=1; i<4; ++i) {
    for (int j=1; j<4; ++j) {
      cout << ' ' << h_Dphi_Dy_jj.bin(i,j).sig;
    }
    cout << endl;
  }

  return 0;
}

