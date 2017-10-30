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
bool is_mc;
// ==================================================================

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
  double bkg = 0, sig = 0; // for significance

  void operator++() noexcept {
    if (is_mc) sig += weight;
    else bkg += weight;
  }
};
double hist_bin::weight;

template <typename... Axes>
using hist = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;

using re_axis = typename re_axes::axis_type;
template <size_t N>
using re_hist = ivanp::binner<hist_bin,
  ivanp::tuple_of_same_t<ivanp::axis_spec<re_axis>,N>>;

// functions applied to variables
inline double phi_pi4(double phi) noexcept {
  phi += M_PI_4;
  return (phi <= M_PI ? phi : phi - M_PI);
}

int main(int argc, const char* argv[]) {
  const std::array<double,2> myy_range{105e3,160e3}, myy_window{121e3,129e3};
  double data_factor = len(myy_window)/(len(myy_range)-len(myy_window));
  double lumi = 0., lumi_in = 0., mc_factor = 1.;

  std::vector<mxaod> mxaods;
  mxaods.reserve(argc-1);
  std::string bins_file("superfine.bins"), fout_name("superfine.root");

  for (int a=1; a<argc; ++a) { // loop over arguments
    // validate args and parse names of input files
    static const std::regex data_re(
      "^(.*/)?data.*_(\\d*)ipb.*\\.root$", std::regex::optimize);
    static const std::regex mc_re(
      "^(.*/)?mc.*\\.root$", std::regex::optimize);
    static const std::regex bins_re(
      "^(.*/)?.*\\.bins$", std::regex::optimize);
    static const std::regex lumi_re(
      "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?) *i([pf])b$",
      std::regex::optimize);
    static const std::regex fout_re(
      "out:(.+\\.root)", std::regex::optimize);
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
    } else if (std::regex_search(arg,end,match,bins_re)) { // MC
      cout << "\033[36mBinning\033[0m: " << arg << endl;
      bins_file = arg;
    } else if (std::regex_search(arg,end,match,fout_re)) { // MC
      fout_name = match[1];
    } else if (std::regex_search(arg,end,match,lumi_re)) {
      lumi = std::stod(match[1]);
      if (arg[match.position(3)]=='f') lumi *= 1e3; // femto to pico
    } else {
      cerr << "arg error: unrecognized argument: " << arg << endl;
      return 1;
    }
  }
  if (!mxaods.size()) {
    cerr << "Must specify at least 1 .root file" << endl;
    return 1;
  }
  cout << "\033[36mOutput file\033[0m: " << fout_name << endl;
  cout << "\n\033[36mTotal data lumi\033[0m: "
       << lumi_in << " ipb" << endl;
  if (lumi!=0.) data_factor *= (lumi / lumi_in);
  else lumi = lumi_in;
  cout << "Scaling to " << lumi << " ipb" << endl << endl;

  // Histogram definitions ==========================================
  re_axes ra(bins_file);
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  h_(pT_yy) h_(yAbs_yy) h_(cosTS_yy) h_(pTt_yy) h_(Dy_y_y)
  h_(HT) h_(HT_yy)
  h_(pT_j1) h_(pT_j2) h_(pT_j3)
  h_(yAbs_j1) h_(yAbs_j2)
  h_(Dphi_j_j) h_(Dphi_j_j_signed)
  h_(Dy_j_j) h_(m_jj)
  h_(pT_yyjj) h_(Dphi_yy_jj)
  h_(sumTau_yyj) h_(maxTau_yyj)
  h_(pT_yy_0j) h_(pT_yy_1j) h_(pT_yy_2j) h_(pT_yy_3j)
  h_(pT_j1_excl)

  h_(xH) h_(x1) h_(x2)

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
        const double n_all = h->GetBinContent(3);
        cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all << endl;
        mc_factor = lumi/n_all;
        break;
      }
    } else { // Data
      hist_bin::weight = data_factor;
    }

    // read variables ===============================================
    TTreeReader reader("CollectionTree",*file);
    optional<TTreeReaderValue<Float_t>> _cs_br_fe, _weight;
    if (is_mc) {
      _cs_br_fe.emplace(reader,"HGamEventInfoAuxDyn.crossSectionBRfilterEff");
      _weight.emplace(reader,"HGamEventInfoAuxDyn.weight");
    }
    TTreeReaderValue<Char_t> _isPassed(reader,"HGamEventInfoAuxDyn.isPassed");

#define VAR_GEN_(NAME, TYPE, STR) \
  TTreeReaderValue<TYPE> _##NAME(reader, "HGamEventInfoAuxDyn." STR);
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
      if (!in(m_yy,myy_range)) continue;

      const bool is_in_window = in(m_yy,myy_window);

      if (is_mc) { // signal from MC
        if (!is_in_window) continue;
        hist_bin::weight = (**_weight) * (**_cs_br_fe) * mc_factor;
      } else { // background from data
        if (is_in_window) continue;
      }

      // FILL HISTOGRAMS ============================================

      const auto nj = *_N_j;

      const auto pT_yy = *_pT_yy*1e-3;
      const auto yAbs_yy = *_yAbs_yy;
      const auto cosTS_yy = std::abs(*_cosTS_yy);
      const auto Dy_y_y = std::abs(*_Dy_y_y);

      h_pT_yy(pT_yy);
      h_yAbs_yy(yAbs_yy);
      h_cosTS_yy(cosTS_yy);

      h_Dy_y_y(Dy_y_y);
      h_pTt_yy(*_pTt_yy*1e-3);

      const auto HT = *_HT*1e-3;
      h_HT(HT);
      h_HT_yy(HT+pT_yy);
      h_xH(pT_yy/HT);

      if (nj == 0) h_pT_yy_0j(pT_yy);

      if (nj < 1) continue; // 1 jet --------------------------------

      const auto pT_j1 = *_pT_j1*1e-3;

      h_pT_j1(pT_j1);

      h_yAbs_j1(*_yAbs_j1);

      h_sumTau_yyj(*_sumTau_yyj*1e-3);
      h_maxTau_yyj(*_maxTau_yyj*1e-3);

      h_x1(pT_j1/HT);

      if (nj == 1) {
        h_pT_j1_excl(pT_j1);
        h_pT_yy_1j(pT_yy);
      }

      if (nj < 2) continue; // 2 jets -------------------------------

      const auto pT_j2   = *_pT_j2*1e-3;
      const auto dphi_jj = std::abs(*_Dphi_j_j);
      const auto   dy_jj = std::abs(*_Dy_j_j);
      const auto    m_jj = *_m_jj*1e-3;

      h_pT_j2(pT_j2);
      h_yAbs_j2(*_yAbs_j2);

      h_Dphi_yy_jj(M_PI - std::abs(*_Dphi_yy_jj));

      h_Dphi_j_j_signed(*_Dphi_j_j_signed);
      h_Dphi_j_j(dphi_jj);
      h_Dy_j_j(dy_jj);
      h_m_jj(m_jj);

      h_pT_yyjj(*_pT_yyjj*1e-3);

      h_x2(pT_j2/HT);

      if (nj == 2) h_pT_yy_2j(pT_yy);

      if (nj < 3) continue; // 3 jets -------------------------------

      h_pT_j3(*_pT_j3*1e-3);
      h_pT_yy_3j(pT_yy);
    }

    file->Close();
  }

  const auto fout = std::make_unique<TFile>(fout_name.c_str(),"recreate");

  for (const auto& h : re_hist<1>::all) {
    const auto& ax = h->axis();
    TH1D *sig = new TH1D((h.name+"_sig").c_str(),"",ax.nbins(),ax.min(),ax.max());
    TH1D *bkg = new TH1D((h.name+"_bkg").c_str(),"",ax.nbins(),ax.min(),ax.max());

    int i = 0;
    for (const auto& bin : h->bins()) {
      sig->SetBinContent(i,bin.sig);
      bkg->SetBinContent(i,bin.bkg);
      ++i;
    }
  }

  fout->Write();

  return 0;
}
