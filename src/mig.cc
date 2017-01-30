#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <memory>
#include <regex>

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

template <typename T>
struct var<TTreeReaderValue<T>> {
  TTreeReaderValue<T> det, truth;
  var(TTreeReader& tr, const std::string& name)
  : det(tr,("HGamEventInfoAuxDyn."+name).c_str()),
    truth(tr,("HGamTruthEventInfoAuxDyn."+name).c_str())
  { }
  using ret_type = std::conditional_t<
    std::is_floating_point<T>::value, double, T>;
  inline var<ret_type> operator*() { return { *det, *truth }; }
  template <typename F>
  inline auto operator()(F f) noexcept(noexcept(f(std::declval<T>())))
  -> var<decltype(f(std::declval<T>()))> { return { f(*det), f(*truth) }; }
};

// global variables =================================================
double n_all_inv, lumi;
var<bool> passed;
// ==================================================================

struct hist_bin {
  static double weight;
  double tmp, x;

  inline hist_bin& operator++() noexcept { tmp += weight; return *this; }

  void merge() noexcept {
    x += tmp*n_all_inv;
    tmp = 0;
  }
};
double hist_bin::weight;

using hist = ivanp::binner<hist_bin, ivanp::tuple_of_same_t<
  ivanp::axis_spec<migration_axis<double>,false,false>,2>>;

auto make_hist(const char* name, const re_axes& ra, unsigned nchecks) {
  migration_axis<double> axis(&*ra[name],nchecks);
  return hist( name, axis, axis );
}

template <typename T, typename... B>
void fill(hist& h, const var<T>& x, const var<B>&... checks) {
  auto hbin = 
  h( std::tie(x.det,   passed.det,   checks.det...),
     std::tie(x.truth, passed.truth, checks.truth...) );
  test( hbin )
}

inline double _abs(double x) noexcept { return std::abs(x); }

int main(int argc, char* argv[]) {
  const std::array<double,2> myy_range{105e3,160e3};
  double lumi = 0.;

  std::vector<std::unique_ptr<TFile>> mxaods;
  mxaods.reserve(argc-2);
  bool lumi_arg = false;
  for (int a=1; a<argc; ++a) { // loop over arguments
    // validate args and parse names of input files
    static const std::regex mc_re(
      "^(.*/)?mc.*\\.root$", std::regex::optimize);
    static const std::regex lumi_re(
      "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?) *i([pf])b$",
      std::regex::optimize);
    std::cmatch match;

    const char *arg = argv[a], *end = arg+std::strlen(arg);
    if (std::regex_search(arg,end,match,mc_re)) { // MC
      cout << "\033[36mMC\033[0m: " << arg << endl;
      mxaods.emplace_back(new TFile(arg,"read"));
      if (mxaods.back()->IsZombie()) return 1;
    } else if (std::regex_search(arg,end,match,lumi_re)) {
      if (lumi_arg) {
        cerr << "arg error: repeated lumi arg: " << arg << endl;
        return 1;
      } else lumi_arg = true;
      lumi = std::stod(match[1]);
      if (arg[match.position(3)]=='f') lumi *= 1e3; // femto to pico
    } else {
      cerr << "arg error: unrecognized argument: " << arg << endl;
      return 1;
    }
  }
  if (!lumi_arg) {
    cerr << "arg error: must specify luminosity" << endl;
    return 1;
  }
  cout << "\n\033[36mScaling to " << lumi << " ipb\033[0m" << endl << endl;

  re_axes ra("hgam.bins");
  // Histogram definitions ==========================================
#define h_(NAME, N) auto h_##NAME = make_hist(#NAME,ra,N);

  // hist<ivanp::index_axis<Int_t>>
  //   h_total("total",{0,1}),
  //   h_N_j_excl("N_j_excl",{0,4}),
  //   h_N_j_incl("N_j_incl",{0,4}),
  //   h_VBF("VBF",{1,4});

  // h_(pT_yy) h_(yAbs_yy) h_(cosTS_yy) h_(pTt_yy) h_(Dy_y_y)
  // h_(HT)
  // h_(pT_j1) h_(pT_j2) h_(pT_j3)
  // h_(yAbs_j1) h_(yAbs_j2)
  // h_(Dphi_j_j) h_(Dphi_j_j_signed)
  // h_(Dy_j_j) h_(m_jj)
  // h_(pT_yyjj) h_(Dphi_yy_jj)
  // h_(sumTau_yyj) h_(maxTau_yyj)
  // h_(pT_yy_0j) h_(pT_yy_1j) h_(pT_yy_2j) h_(pT_yy_3j)
  // h_(pT_j1_excl)

  h_(Dphi_j_j,2)

  for (auto& file : mxaods) { // loop over MxAODs
    cout << "\033[36mMC\033[0m: " << file->GetName() << endl;

    { TIter next(file->GetListOfKeys());
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

    // read variables ===============================================
    TTreeReader reader("CollectionTree",file.get());
    TTreeReaderValue<Char_t> _isPassed(reader,"HGamEventInfoAuxDyn.isPassed");
    TTreeReaderValue<Float_t> _cs_br_fe(reader,
      "HGamEventInfoAuxDyn.crossSectionBRfilterEff");
    TTreeReaderValue<Float_t> _weight(reader,
      "HGamEventInfoAuxDyn.weight");
    TTreeReaderValue<Char_t> _isFiducial(reader,
      "HGamTruthEventInfoAuxDyn.isFiducial");

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

      // diphoton mass cut
      const auto m_yy = *_m_yy;
      if (!in(m_yy.det,myy_range)) continue;

      hist_bin::weight = (*_weight)*(*_cs_br_fe);
      test( hist_bin::weight )

      passed = {*_isPassed,*_isFiducial && in(m_yy.truth,myy_range)};

      // FILL HISTOGRAMS ============================================
      const auto nj = *_N_j;

      const auto dphi_jj = _Dphi_j_j(_abs);

      if (dphi_jj.det == 99) continue;

      test( nj.det )
      test( nj.truth )
      test( dphi_jj.det )
      test( dphi_jj.truth )

      fill(h_Dphi_j_j, dphi_jj, nj>=2);

      // if (ent==4) break;
      // cout << endl;
    }

    for (const auto& h : hist::all)
      for (auto& b : h->bins()) b.merge();

    file->Close();
  }

  for (const auto& h : hist::all) {
    cout << h.name << endl;
    for (auto& b : h->bins()) cout << b.x*lumi << endl;
    cout << endl;
  }

  return 0;
}
