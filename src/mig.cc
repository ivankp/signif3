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

#define VAR_ALWAYS_MC
#include "truth_reco_var.hh"

// global variables =================================================
double n_all_inv, lumi;
var<bool> passed;
// ==================================================================

struct hist_bin {
  static double weight;
  double tmp, w;

  inline hist_bin& operator++() noexcept { tmp += weight; return *this; }

  void merge() noexcept {
    w += tmp*n_all_inv;
    tmp = 0;
  }
};
double hist_bin::weight;

template <typename T>
using hist = ivanp::binner<hist_bin, ivanp::tuple_of_same_t<
  ivanp::axis_spec<migration_axis<T>,true,true>,2>>;
  // ivanp::axis_spec<migration_axis<T>,false,false>,2>>;

auto make_hist(const char* name, const re_axes& ra, unsigned nchecks) {
  migration_axis<double> axis(&*ra[name],nchecks);
  return hist<double>( name, axis, axis );
}

template <typename A, typename T, typename... B>
void fill(ivanp::binner<hist_bin,std::tuple<A,A>>& h,
  const var<T>& x, const var<B>&... checks
) {
  h( std::tie(x.det,   passed.det,   checks.det...),
     std::tie(x.truth, passed.truth, checks.truth...) );
}

template <typename A>
std::ostream& operator<<(std::ostream& o,
  const ivanp::named_ptr<ivanp::binner<hist_bin,std::tuple<A,A>>>& h
) {
  std::ios::fmtflags f(cout.flags());
  auto prec = cout.precision();
  cout.precision(2);
  cout << std::fixed;

  o << "\033[32m" << h.name << "\033[0m\n";
  const unsigned n = h->axis().nbins() + A::nover::value;
  unsigned i = 0;
  for (auto& b : h->bins()) {
    cout << ' ' << std::setw(7) << b.w*lumi;
    if (++i == n) { cout << endl; i = 0; }
  }

  cout.flags(f);
  cout.precision(prec);
  return o;
}

int main(int argc, char* argv[]) {
  if (argc==1) {
    cout << "usage: " << argv[0] << " mc*.root ?i[pf]b" << endl;
    return 1;
  }
  const std::array<double,2> myy_range{105e3,160e3};
  // const std::array<double,2> myy_range{121e3,129e3};

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

  ivanp::index_axis<Int_t,true> nj_axis(0,4);
  hist<Int_t>
    h_N_j_excl("N_j_excl",{&nj_axis,1},{&nj_axis,1});

  h_(pT_yy,1) h_(yAbs_yy,1) h_(cosTS_yy,1) h_(pTt_yy,1) h_(Dy_y_y,1)
  h_(HT,1)
  h_(pT_j1,2) h_(pT_j2,2) h_(pT_j3,2)
  h_(yAbs_j1,2) h_(yAbs_j2,2)
  h_(Dphi_j_j,2) h_(Dphi_j_j_signed,2)
  h_(Dy_j_j,2) h_(m_jj,2)
  h_(pT_yyjj,2) h_(Dphi_yy_jj,2)
  h_(sumTau_yyj,2) h_(maxTau_yyj,2)
  h_(pT_yy_0j,2) h_(pT_yy_1j,2) h_(pT_yy_2j,2) h_(pT_yy_3j,2)
  h_(pT_j1_excl,2)

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
      // test( hist_bin::weight )

      passed = {*_isPassed,*_isFiducial && in(m_yy.truth,myy_range)};

      // FILL HISTOGRAMS ============================================
      const auto nj = *_N_j;

      const auto pT_yy = _pT_yy/1e3;
      fill(h_pT_yy, pT_yy);
      fill(h_pT_yy_0j, pT_yy, nj==0);
      fill(h_pT_yy_1j, pT_yy, nj==1);
      fill(h_pT_yy_2j, pT_yy, nj==2);
      fill(h_pT_yy_3j, pT_yy, nj>=3);

      fill(h_yAbs_yy, *_yAbs_yy);
      fill(h_cosTS_yy, abs(_cosTS_yy));

      fill(h_Dy_y_y, abs(_Dy_y_y));
      fill(h_pTt_yy, _pTt_yy/1e3);

      fill(h_N_j_excl, nj);

      fill(h_HT, _HT/1e3);

      const auto pT_j1 = _pT_j1/1e3;
      fill(h_pT_j1, pT_j1, nj>=1);
      fill(h_pT_j1_excl, pT_j1, nj==1);
      fill(h_yAbs_j1, *_yAbs_j1, nj>=1);

      fill(h_pT_j2, _pT_j2/1e3, nj>=2);
      fill(h_yAbs_j2, *_yAbs_j2, nj>=2);

      fill(h_pT_j3, _pT_j3/1e3, nj>=3);

      fill(h_sumTau_yyj, _sumTau_yyj/1e3, nj>=1);
      fill(h_maxTau_yyj, _maxTau_yyj/1e3, nj>=1);

      fill(h_Dphi_j_j, abs(_Dphi_j_j), nj>=2);
      fill(h_Dy_j_j, abs(_Dy_j_j), nj>=2);

      fill(h_Dphi_yy_jj, _Dphi_yy_jj%[](auto x){ return M_PI - std::abs(x);},
           nj>=2);

      fill(h_Dphi_j_j_signed, *_Dphi_j_j_signed, nj>=2);
      fill(h_Dphi_j_j, abs(_Dphi_j_j), nj>=2);
      fill(h_m_jj, _m_jj/1e3, nj>=2);

      fill(h_pT_yyjj, _pT_yyjj/1e3, nj>=2);
    }

    for (const auto& h : hist<Int_t>::all)
      for (auto& b : h->bins()) b.merge();
    for (const auto& h : hist<double>::all)
      for (auto& b : h->bins()) b.merge();

    file->Close();
  }

  cout << endl;
  for (const auto& h : hist<Int_t>::all) cout << h << endl;
  for (const auto& h : hist<double>::all) cout << h << endl;

  return 0;
}
