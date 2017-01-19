#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <memory>
#include <regex>
#include <experimental/optional>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TKey.h>

#include "binner.hh"
#include "re_axes.hh"
#include "catstr.hh"
#include "timed_counter.hh"
#include "array_ops.hh"
#include "exception.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

bool is_mc, is_fiducial; // global variables
double n_all_inv, fw, lumi;

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
      bkg += tmp*n_all_inv;
    }
    tmp = 0;
  }
  void compute_signif() noexcept {
    bkg *= fw;
    sig *= lumi;
    const double sb = sig + bkg;
    signif = sig/std::sqrt(sb);
    sb_purity = sig/sb;
    tr_purity = sig_truth/sig;
  };
};
double hist_bin::weight;

template <typename... Axes>
using hist = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;

using re_axis = typename re_axes::axis_type;
template <size_t N>
using re_hist = ivanp::binner<hist_bin,
  ivanp::tuple_of_same_t<ivanp::axis_spec<re_axis>,N>>;

template <typename... Axes>
std::ostream& operator<<(std::ostream& o,
  const ivanp::named_ptr<hist<Axes...>>& h
) {
  const auto prec = o.precision();
  const std::ios::fmtflags f( o.flags() );
  o << "\033[32m" << h.name << "\033[0m\n";
  const auto& a = h->axis();
  for (unsigned i=1, n=a.nbins()+2; i<n; ++i) {
    const auto& b = h->bin(i);
    o << "\033[35m[" << a.lower(i) << ',';
    if (i==n-1) o << "∞";
    else o << a.upper(i);
    o << ")\033[0m "
      << round(b.sig) << ' '
      << round(b.bkg) << ' '
      << std::setprecision(2) << std::fixed
      << b.signif << ' '
      << 100*b.sb_purity << "% "
      << 100*b.tr_purity << '%'
      << std::setprecision(prec) << endl;
    o.flags( f );
  }
  return o;
}

template <typename T>
struct var {
  T det, truth;
  inline bool operator==(const T& x) const noexcept { return (det == x); }
  inline bool operator< (const T& x) const noexcept { return (det <  x); }
  inline bool operator> (const T& x) const noexcept { return (det >  x); }
};
template <typename T>
struct var<TTreeReaderValue<T>> {
  TTreeReaderValue<T> det;
  std::experimental::optional<TTreeReaderValue<T>> truth;
  var(TTreeReader& tr, const std::string& name)
  : det(tr,("HGamEventInfoAuxDyn."+name).c_str())
  {
    if (is_mc) truth.emplace(tr,("HGamTruthEventInfoAuxDyn."+name).c_str());
  }
  using ret_type = std::conditional_t<
    std::is_floating_point<T>::value, double, T>;
  var<ret_type> operator*() {
    if (truth) return { *det, **truth };
    else return { *det, { } };
  }
  template <typename F>
  var<ret_type> operator()(F f) {
    if (truth) return { f(*det), f(**truth) };
    else return { f(*det), { } };
  }
};

template <typename T, typename... Axes>
void fill(hist<Axes...>& h, const var<T>& x) {
  if (is_mc) {
    const auto bin_det = h.find_bin(x.det);
    const auto bin_truth = h.find_bin(x.truth);
    h.fill_bin(bin_det, is_fiducial && (bin_det == bin_truth));
  } else h(x.det);
}

double d1e3(double x) { return x*1e-3; }
double _abs(double x) { return std::abs(x); }

int main(int argc, const char* argv[])
{
  const std::array<double,2> mass_range{105e3,160e3}, mass_window{121e3,129e3};
  fw = len(mass_window)/(len(mass_range)-len(mass_window));
  lumi = 0.;

  re_axes ra("hgam.bins");
// #define a_(name) auto a_##name = ra[#name];
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  hist<ivanp::index_axis<Int_t>>
    h_total("total",{0,1}),
    h_N_j_excl("N_j_excl",{0,4}),
    h_N_j_incl("N_j_incl",{0,4}),
    h_VBF("VBF",{0,3});

  h_(pT_yy) h_(yAbs_yy) h_(cosTS_yy) h_(pTt_yy) h_(Dy_y_y)
  h_(HT)
  h_(pT_j1) h_(pT_j2) h_(pT_j3)
  h_(yAbs_j1) h_(yAbs_j2)
  h_(Dphi_j_j) h_(Dphi_j_j_signed)
  h_(Dy_j_j) h_(m_jj)
  h_(pT_yyjj) h_(Dphi_yy_jj)
  h_(sumTau_yyj) h_(maxTau_yyj)
  h_(pT_yy_0j) h_(pT_yy_1j) h_(pT_yy_2j) h_(pT_yy_3j)
  h_(pT_j1_excl)

  // h_Dphi_Dy_jj("Dphi_Dy_jj",{0.,M_PI_2,M_PI},{0.,2.,8.8}),
  // h_Dphi_pi4_Dy_jj("Dphi_pi4_Dy_jj",{0.,M_PI_2,M_PI},{0.,2.,8.8}),
  // h_cosTS_pT_yy("cosTS_pT_yy",{0.,0.5,1.},{0.,30.,120.,400.}),
  // h_pT_yy_pT_j1("pT_yy_pT_j1",{0.,30.,120.,400.},{30.,65.,400.});

  for (int f=1; f<argc; ++f) {
    { static const std::regex data_re(".*/data.*_(\\d*)ipb.*\\.root$");
      static const std::regex mc_re(".*/mc.*\\.root$");
      std::cmatch match;

      const char *fname = argv[f], *end = fname+std::strlen(fname);
      if (std::regex_search(fname,end,match,data_re)) {
        is_mc = false;
        hist_bin::weight = 1;
        const auto flumi = std::stod(match[1]);
        lumi += flumi;
        cout << "\033[36mData\033[0m: " << fname << endl;
        cout << "\033[36mLumi\033[0m: " << flumi << " ipb" << endl;
      } else if (std::regex_search(fname,end,match,mc_re)) {
        is_mc = true;
        cout << "\033[36mMC\033[0m: " << fname << endl;
      } else throw ivanp::exception("Unexpected file name: ",fname);
    }

    auto file = std::make_unique<TFile>(argv[f],"read");
    if (file->IsZombie()) return 1;

    n_all_inv = 1;
    if (is_mc) {
      TIter next(file->GetListOfKeys());
      TKey *key;
      while ((key = static_cast<TKey*>(next()))) {
        std::string name(key->GetName());
        if (name.substr(0,8)!="CutFlow_" ||
            name.substr(name.size()-18)!="_noDalitz_weighted") continue;
        TH1 *h = static_cast<TH1*>(key->ReadObj());
        cout << h->GetName() << endl;
        n_all_inv = h->GetBinContent(3);
        cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all_inv << endl;
        n_all_inv = 1./n_all_inv;
        break;
      }
    }

    TTreeReader reader("CollectionTree",file.get());
    std::experimental::optional<TTreeReaderValue<Float_t>> _cs_br_fe, _weight;
    if (is_mc) {
      _cs_br_fe.emplace(reader,"HGamEventInfoAuxDyn.crossSectionBRfilterEff");
      _weight .emplace(reader,"HGamEventInfoAuxDyn.weight");
    }
    TTreeReaderValue<Char_t> _isPassed(reader,"HGamEventInfoAuxDyn.isPassed");
    TTreeReaderValue<Char_t> _isFiducial(reader,"HGamTruthEventInfoAuxDyn.isFiducial");

#define VAR_GEN_(NAME, TYPE, STR) \
  var<TTreeReaderValue<TYPE>> _##NAME(reader, STR);
#define VAR_(NAME) VAR_GEN_(NAME, Float_t, #NAME)
#define VAR30_(NAME) VAR_GEN_(NAME, Float_t, #NAME "_30")

    VAR_GEN_(N_j, Int_t, "N_j_30")

    VAR_(m_yy) VAR_(pT_yy) VAR_(yAbs_yy) VAR_(cosTS_yy) VAR_(pTt_yy)
    VAR_(Dy_y_y)

    VAR30_(HT)
    VAR30_(pT_j1)      VAR30_(pT_j2)      VAR30_(pT_j3)
    VAR30_(yAbs_j1)    VAR30_(yAbs_j2)
    VAR30_(Dphi_j_j)   VAR_GEN_(Dphi_j_j_signed,Float_t,"Dphi_j_j_30_signed")
    VAR30_(Dy_j_j)     VAR30_(m_jj)
    VAR30_(sumTau_yyj) VAR30_(maxTau_yyj)
    VAR30_(pT_yyjj)    VAR30_(Dphi_yy_jj)

    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
      if (!*_isPassed) continue;

      const auto m_yy = *_m_yy;
      if (!in(m_yy.det,mass_range)) continue;

      if (is_mc) { // signal from MC
        if (!in(m_yy.det,mass_window)) continue;
        hist_bin::weight = (**_weight)*(**_cs_br_fe);
      } else { // background from data
        if (in(m_yy.det,mass_window)) continue;
      }

      is_fiducial = *_isFiducial;

      // FILL HISTOGRAMS =====================================================

      const auto nj = *_N_j;

      const auto pT_yy = _pT_yy(d1e3);
      const auto yAbs_yy = *_yAbs_yy;
      const auto cosTS_yy = _cosTS_yy(_abs);
      const auto Dy_y_y = _Dy_y_y(_abs);

      h_total(0, is_fiducial);

      fill(h_pT_yy, pT_yy);
      fill(h_yAbs_yy, yAbs_yy);
      fill(h_cosTS_yy, cosTS_yy);

      fill(h_Dy_y_y, Dy_y_y);
      fill(h_pTt_yy, _pTt_yy(d1e3));
      // h_cosTS_pT_yy( cosTS_yy, pT_yy );

      fill(h_N_j_excl, nj);

      fill(h_HT, _HT(d1e3));

      if (nj == 0) fill(h_pT_yy_0j, pT_yy);

      if (nj < 1) continue; // 1 jet -------------------------------

      const auto pT_j1 = _pT_j1(d1e3);

      fill(h_pT_j1, pT_j1);

      fill(h_yAbs_j1, *_yAbs_j1);

      fill(h_sumTau_yyj, _sumTau_yyj(d1e3));
      fill(h_maxTau_yyj, _maxTau_yyj(d1e3));

      if (nj == 1) {
        fill(h_pT_j1_excl, pT_j1);
        fill(h_pT_yy_1j, pT_yy);
      }

      // h_pT_yy_pT_j1( pT_yy, pT_j1 );

      if (nj < 2) continue; // 2 jets ------------------------------

      const auto dphi_jj = _Dphi_j_j(_abs);
      const auto   dy_jj = _Dy_j_j(_abs);
      const auto    m_jj = _m_jj(d1e3);

      fill(h_pT_j2, _pT_j2(d1e3));
      fill(h_yAbs_j2, *_yAbs_j2);

      fill(h_Dphi_yy_jj, _Dphi_yy_jj([](auto x){ return M_PI - std::abs(x);}));

      fill(h_Dphi_j_j_signed, *_Dphi_j_j_signed);
      fill(h_Dphi_j_j, dphi_jj);
      fill(h_Dy_j_j, dy_jj);
      fill(h_m_jj, m_jj);

      fill(h_pT_yyjj, _pT_yyjj(d1e3));

      if (nj == 2) fill(h_pT_yy_2j, pT_yy);

      // VBF --------------------------------------------------------
      var<double> pT_j3{0.,0.};
      if (nj > 2) pT_j3 = _pT_j3(d1e3);

      if ( (m_jj > 600.) && (dy_jj > 4.0) ) {
        if (pT_j3 < 30.) h_VBF(0);
        if (pT_j3 < 25.) h_VBF(1);
      }
      if ( (m_jj > 400.) && (dy_jj > 2.8) ) {
        if (pT_j3 < 30.) h_VBF(2);
      }
      // ------------------------------------------------------------

      if (nj < 3) continue; // 3 jets ------------------------------

      fill(h_pT_yy_3j, pT_yy);
      fill(h_pT_j3, pT_j3);

      // h_Dphi_Dy_jj( dphi_jj, dy_jj );
      // h_Dphi_pi4_Dy_jj( phi_pi4(dphi_jj), dy_jj );
    }

    for (const auto& h : hist<ivanp::index_axis<Int_t>>::all)
      for (auto& b : h->bins()) b.merge();
    for (const auto& h : re_hist<1>::all)
      for (auto& b : h->bins()) b.merge();

  }

  h_N_j_incl = h_N_j_excl;
  h_N_j_incl.integrate_left();

  cout << "\n\033[36mTotal lumi\033[0m: " << lumi <<" ipb\n"<< endl;

  for (const auto& h : hist<ivanp::index_axis<Int_t>>::all) {
    for (auto& b : h->bins()) b.compute_signif();
    cout << h << endl;
  }
  for (const auto& h : re_hist<1>::all) {
    for (auto& b : h->bins()) b.compute_signif();
    cout << h << endl;
  }

  return 0;
}
