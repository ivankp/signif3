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
using std::experimental::optional;

// global variables =================================================
bool is_mc, is_fiducial, is_in_window;
double n_all_inv, factor, lumi = 0.;
// ==================================================================

template <typename T>
struct var {
  T det, truth;
  inline bool operator==(const T& x) const noexcept { return (det == x); }
  inline bool operator< (const T& x) const noexcept { return (det <  x); }
  inline bool operator> (const T& x) const noexcept { return (det >  x); }

  template <typename F>
  inline auto operator()(F f) const noexcept(noexcept(f(std::declval<T>())))
  -> var<decltype(f(std::declval<T>()))> {
    if (is_mc) return { f(det), f(truth) };
    else return { f(det), { } };
  }
};

template <typename T>
struct var<TTreeReaderValue<T>> {
  TTreeReaderValue<T> det;
  optional<TTreeReaderValue<T>> truth;
  var(TTreeReader& tr, const std::string& name)
  : det(tr,("HGamEventInfoAuxDyn."+name).c_str())
  {
    if (is_mc) truth.emplace(tr,("HGamTruthEventInfoAuxDyn."+name).c_str());
  }
  using ret_type = std::conditional_t<
    std::is_floating_point<T>::value, double, T>;
  inline var<ret_type> operator*() {
    if (truth) return { *det, **truth };
    else return { *det, { } };
  }
  template <typename F>
  inline auto operator()(F f) noexcept(noexcept(f(std::declval<T>())))
  -> var<decltype(f(std::declval<T>()))> {
    if (truth) return { f(*det), f(**truth) };
    else return { f(*det), { } };
  }
};

class mxaod {
  TFile *ptr;
  bool _is_mc;
public:
  mxaod(const char* fname, bool _is_mc)
  : ptr(new TFile(fname,"read")), _is_mc(_is_mc) {
    if (ptr->IsZombie()) throw ivanp::exception("cannot open file ",fname);
  }
  ~mxaod() { delete ptr; }
  mxaod(const mxaod& o) = delete;
  mxaod(mxaod&& o) noexcept: ptr(o.ptr), _is_mc(o._is_mc) { o.ptr = nullptr; }
  inline TFile* operator->() const noexcept { return ptr; }
  inline TFile* operator* () const noexcept { return ptr; }
  inline bool is_mc() const noexcept { return _is_mc; }
};

struct hist_bin {
  static double weight;
  double
    tmp_sig, // temporary for significance
             // in window for mc, in range but not window for data
    tmp_reco,  // mc reco in range, for purity
    tmp_truth, // mc truth in range, for purity
    bkg, sig,    // for significance
    reco, truth; // for purity

  hist_bin(): tmp_sig(0), tmp_reco(0), tmp_truth(0),
              bkg(0), sig(0), reco(0), truth(0) { }

  inline void operator()(bool truth_match=true) noexcept {
    if (is_mc) {
      if (is_in_window) tmp_sig += weight; // cut for significance
      tmp_reco += weight;
      // is_fiducial includes mass check
      if (is_fiducial && truth_match) tmp_truth += weight;
    } else {
      // alway fill data here
      // the cut is in the event loop
      tmp_sig += weight;
    }
  }

  void merge() noexcept {
    if (is_mc) {
      sig += tmp_sig*n_all_inv; tmp_sig = 0;
      reco += tmp_reco*n_all_inv; tmp_reco = 0;
      truth += tmp_truth*n_all_inv; tmp_truth = 0;
    } else {
      bkg += tmp_sig; tmp_sig = 0;
    }
  }
};
double hist_bin::weight;

std::ostream& operator<<(std::ostream& o, const hist_bin& b) {
  const double // compute significance and purity
    bkg = b.bkg * factor,
    sig = b.sig * lumi,
    reco = b.reco * lumi,
    truth = b.truth * lumi,
    signif = sig/std::sqrt(sig+bkg),
    purity = truth/reco;

  const auto prec = o.precision();
  const std::ios::fmtflags f( o.flags() );
  o << round(sig) << ' '
    << round(bkg) << ' '
    << std::setprecision(2) << std::fixed
    << signif << ' '
    << (100*sig/(sig+bkg)) << "% "
    << (100*purity) << '%'
    << std::setprecision(prec);
  o.flags( f );
  return o;
}

template <typename... Axes>
using hist = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;

using re_axis = typename re_axes::axis_type;
template <size_t N>
using re_hist = ivanp::binner<hist_bin,
  ivanp::tuple_of_same_t<ivanp::axis_spec<re_axis>,N>>;

template <typename A1>
std::ostream& operator<<(std::ostream& o,
  const ivanp::named_ptr<hist<A1>>& h
) {
  o << "\033[32m" << h.name << "\033[0m\n";
  const auto& a = h->axis();
  for (unsigned i=1, n=a.nbins()+2; i<n; ++i) {
    o << "\033[35m[" << a.lower(i) << ',';
    if (i==n-1) o << "∞";
    else o << a.upper(i);
    o << ")\033[0m " << h->bin(i) << endl;
  }
  return o;
}

template <typename A1, typename A2>
std::ostream& operator<<(std::ostream& o,
  const ivanp::named_ptr<hist<A1,A2>>& h
) {
  o << "\033[32m" << h.name << "\033[0m\n";
  const auto& a1 = std::get<0>(h->axes());
  const auto& a2 = std::get<1>(h->axes());
  for (unsigned i1=1, n1=a1.nbins()+2; i1<n1; ++i1) {
    o << "\033[35m[" << a1.lower(i1) << ',';
    if (i1==n1-1) o << "∞";
    else o << a1.upper(i1);
    o << ")\033[0m \n";
    for (unsigned i2=1, n2=a2.nbins()+2; i2<n2; ++i2) {
      o << "  \033[35m[" << a2.lower(i2) << ',';
      if (i2==n2-1) o << "∞";
      else o << a2.upper(i2);
      o << ")\033[0m " << h->bin(i1,i2) << endl;
    }
  }
  return o;
}

template <typename T, typename Axis>
void fill(hist<Axis>& h, const var<T>& x, bool extra_truth_match=true) {
  const auto bin_det = h.find_bin(x.det);
  if (is_mc) {
    const auto bin_truth = h.find_bin(x.truth);
    h.fill_bin(bin_det, (bin_det == bin_truth) && extra_truth_match);
  } else h.fill_bin(bin_det);
}

template <typename T, typename Axis>
void fill_incl(hist<Axis>& h, const var<T>& x) {
  const auto bin_det = h.find_bin(x.det);
  if (is_mc) {
    const auto bin_truth = h.find_bin(x.truth);
    for (unsigned i=bin_det; i!=0; --i)
      h.fill_bin(i, bin_truth >= i);
  } else for (unsigned i=bin_det; i!=0; --i) h.fill_bin(i);
}

template <typename T1, typename T2, typename A1, typename A2>
void fill(hist<A1,A2>& h, const var<T1>& x1, const var<T2>& x2,
  bool extra_truth_match=true
) {
  const auto bin_det = h.find_bin(x1.det,x2.det);
  if (is_mc) {
    const auto bin_truth = h.find_bin(x1.truth,x2.truth);
    h.fill_bin(bin_det, (bin_det == bin_truth) && extra_truth_match);
  } else h.fill_bin(bin_det);
}

// functions applied to variables
inline double d1e3(double x) noexcept { return x*1e-3; }
inline double _abs(double x) noexcept { return std::abs(x); }
inline double phi_pi4(double phi) noexcept {
  phi += M_PI_4;
  return (phi <= M_PI ? phi : phi - M_PI);
}

template <typename F, typename... T>
auto apply(F f, const var<T>&... vars) -> var<decltype(f(vars.det...))> {
  if (is_mc) return { f(vars.det...), f(vars.truth...) };
  else return { f(vars.det...), { } };
}

int main(int argc, const char* argv[]) {
  const std::array<double,2> mass_range{105e3,160e3}, mass_window{121e3,129e3};
  factor = len(mass_window)/(len(mass_range)-len(mass_window));
  double lumi_in = 0.;

  re_axes ra("hgam.bins");
  // Histogram definitions ==========================================
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  hist<ivanp::index_axis<Int_t>>
    h_total("total",{0,1}),
    h_N_j_excl("N_j_excl",{0,4}),
    h_N_j_incl("N_j_incl",{0,4}),
    h_VBF("VBF",{1,4});

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

  using hist2 = hist<
    ivanp::container_axis<std::vector<double>>,
    ivanp::container_axis<std::vector<double>> >;
  hist2
    h_Dphi_Dy_jj("Dphi_Dy_jj",{0.,M_PI_2,M_PI},{0.,2.,8.8}),
    h_Dphi_pi4_Dy_jj("Dphi_pi4_Dy_jj",{0.,M_PI_2,M_PI},{0.,2.,8.8}),
    h_cosTS_pT_yy("cosTS_pT_yy",{0.,0.5,1.},{0.,30.,120.,400.}),
    h_pT_yy_pT_j1("pT_yy_pT_j1",{0.,30.,120.,400.},{30.,65.,400.});

  std::vector<mxaod> mxaods;
  mxaods.reserve(argc-1);
  for (int a=1; a<argc; ++a) { // loop over arguments
    // validate args and parse names of input files
    static const std::regex data_re(
      "^(.*/)?data.*_(\\d*)ipb.*\\.root$", std::regex::optimize);
    static const std::regex mc_re(
      "^(.*/)?mc.*\\.root$", std::regex::optimize);
    static const std::regex lumi_re(
      "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?) *i([pf])b$",
      std::regex::optimize);
    std::cmatch match;

    const char *arg = argv[a], *end = arg+std::strlen(arg);
    if (std::regex_search(arg,end,match,data_re)) { // Data
      const double flumi = std::stod(match[2]);
      lumi_in += flumi;
      cout << "\033[36mData\033[0m: " << arg << endl;
      cout << "\033[36mLumi\033[0m: " << flumi << " ipb" << endl;
      mxaods.emplace_back(arg,false);
    } else if (std::regex_search(arg,end,match,mc_re)) { // MC
      cout << "\033[36mMC\033[0m: " << arg << endl;
      mxaods.emplace_back(arg,true);
    } else if (std::regex_search(arg,end,match,lumi_re)) {
      lumi = std::stod(match[1]);
      if (arg[match.position(3)]=='f') lumi *= 1e3; // femto to pico
    } else throw ivanp::exception("unrecognized argument: ",arg);
  }
  cout << "\n\033[36mTotal data lumi\033[0m: "
       << lumi_in << " ipb" << endl;
  if (lumi!=0.) factor *= (lumi / lumi_in);
  else lumi = lumi_in;
  cout << "Scaling to " << lumi << " ipb" << endl << endl;

  for (auto& file : mxaods) { // loop over MxAODs
    is_mc = file.is_mc();
    cout << "\033[36m" << (is_mc ? "MC" : "Data") << "\033[0m: "
         << file->GetName() << endl;

    if (is_mc) { // MC
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
    } else { // Data
      hist_bin::weight = 1;
      n_all_inv = 1;
    }

    // read variables ===============================================
    TTreeReader reader("CollectionTree",*file);
    optional<TTreeReaderValue<Float_t>> _cs_br_fe, _weight;
    optional<TTreeReaderValue<Char_t>> _isFiducial;
    if (is_mc) {
      _cs_br_fe.emplace(reader,"HGamEventInfoAuxDyn.crossSectionBRfilterEff");
      _weight.emplace(reader,"HGamEventInfoAuxDyn.weight");
      _isFiducial.emplace(reader,"HGamTruthEventInfoAuxDyn.isFiducial");
    }
    TTreeReaderValue<Char_t> _isPassed(reader,"HGamEventInfoAuxDyn.isPassed");

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

    // loop over events =============================================
    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {

      // selection cut
      if (!*_isPassed) continue;

      // diphoton mass cut
      const auto m_yy = *_m_yy;
      if (!in(m_yy.det,mass_range)) continue;

      is_in_window = in(m_yy.det,mass_window);

      if (is_mc) { // signal from MC
        hist_bin::weight = (**_weight)*(**_cs_br_fe);
        is_fiducial = **_isFiducial && in(m_yy.truth,mass_range);
      } else { // background from data
        if (is_in_window) continue;
      }

      // FILL HISTOGRAMS ============================================

      const auto nj = *_N_j;
      bool match_truth_nj;

      const auto pT_yy = _pT_yy(d1e3);
      const auto yAbs_yy = *_yAbs_yy;
      const auto cosTS_yy = _cosTS_yy(_abs);
      const auto Dy_y_y = _Dy_y_y(_abs);

      h_total(0);

      fill(h_pT_yy, pT_yy);
      fill(h_yAbs_yy, yAbs_yy);
      fill(h_cosTS_yy, cosTS_yy);

      fill(h_Dy_y_y, Dy_y_y);
      fill(h_pTt_yy, _pTt_yy(d1e3));
      fill(h_cosTS_pT_yy, cosTS_yy, pT_yy);

      fill(h_N_j_excl, nj);
      fill_incl(h_N_j_incl, nj);

      fill(h_HT, _HT(d1e3));

      if (nj == 0) fill(h_pT_yy_0j, pT_yy, nj.truth==0);

      if (nj < 1) continue; // 1 jet --------------------------------

      match_truth_nj = nj.truth>=1;

      const auto pT_j1 = _pT_j1(d1e3);

      fill(h_pT_j1, pT_j1, match_truth_nj);

      fill(h_yAbs_j1, *_yAbs_j1, match_truth_nj);

      fill(h_sumTau_yyj, _sumTau_yyj(d1e3), match_truth_nj);
      fill(h_maxTau_yyj, _maxTau_yyj(d1e3), match_truth_nj);

      fill(h_pT_yy_pT_j1, pT_yy, pT_j1, match_truth_nj);

      if (nj == 1) {
        match_truth_nj = nj.truth==1;
        fill(h_pT_j1_excl, pT_j1, match_truth_nj);
        fill(h_pT_yy_1j, pT_yy, match_truth_nj);
      }

      if (nj < 2) continue; // 2 jets -------------------------------

      match_truth_nj = nj.truth>=2;

      const auto dphi_jj = _Dphi_j_j(_abs);
      const auto   dy_jj = _Dy_j_j(_abs);
      const auto    m_jj = _m_jj(d1e3);

      fill(h_pT_j2, _pT_j2(d1e3), match_truth_nj);
      fill(h_yAbs_j2, *_yAbs_j2, match_truth_nj);

      fill(h_Dphi_yy_jj, _Dphi_yy_jj([](auto x){ return M_PI - std::abs(x);}),
        match_truth_nj);

      fill(h_Dphi_j_j_signed, *_Dphi_j_j_signed, match_truth_nj);
      fill(h_Dphi_j_j, dphi_jj, match_truth_nj);
      fill(h_Dy_j_j, dy_jj, match_truth_nj);
      fill(h_m_jj, m_jj, match_truth_nj);

      fill(h_pT_yyjj, _pT_yyjj(d1e3), match_truth_nj);

      fill(h_Dphi_Dy_jj, dphi_jj, dy_jj, match_truth_nj);
      fill(h_Dphi_pi4_Dy_jj, dphi_jj(phi_pi4), dy_jj, match_truth_nj);

      if (nj == 2) fill(h_pT_yy_2j, pT_yy, nj.truth==2);

      // VBF --------------------------------------------------------
      var<double> pT_j3{0.,0.};
      if (nj > 2) pT_j3 = _pT_j3(d1e3);

      auto VBF1 = apply([](double m_jj, double dy_jj, double pT_j3) {
        return (m_jj > 600.) && (dy_jj > 4.0) && (pT_j3 < 30.);
      }, m_jj, dy_jj, pT_j3);
      auto VBF2 = apply([](double m_jj, double dy_jj, double pT_j3) {
        return (m_jj > 600.) && (dy_jj > 4.0) && (pT_j3 < 25.);
      }, m_jj, dy_jj, pT_j3);
      auto VBF3 = apply([](double m_jj, double dy_jj, double pT_j3) {
        return (m_jj > 400.) && (dy_jj > 2.8) && (pT_j3 < 30.);
      }, m_jj, dy_jj, pT_j3);

      if (VBF1.det) h_VBF.fill_bin(1,VBF1.det==VBF1.truth);
      if (VBF2.det) h_VBF.fill_bin(2,VBF2.det==VBF2.truth);
      if (VBF3.det) h_VBF.fill_bin(3,VBF3.det==VBF3.truth);
      // ------------------------------------------------------------

      if (nj < 3) continue; // 3 jets -------------------------------

      match_truth_nj = nj.truth>=3;

      fill(h_pT_yy_3j, pT_yy, match_truth_nj);
      fill(h_pT_j3, pT_j3, match_truth_nj);
    }

    for (const auto& h : hist<ivanp::index_axis<Int_t>>::all)
      for (auto& b : h->bins()) b.merge();
    for (const auto& h : re_hist<1>::all)
      for (auto& b : h->bins()) b.merge();
    for (const auto& h : hist2::all)
      for (auto& b : h->bins()) b.merge();

    file->Close();
  }

  for (const auto& h : hist<ivanp::index_axis<Int_t>>::all) cout << h << endl;
  for (const auto& h : re_hist<1>::all) cout << h << endl;
  for (const auto& h : hist2::all) cout << h << endl;

  return 0;
}
