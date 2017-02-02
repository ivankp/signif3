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
bool is_mc;
double n_all_inv, lumi = 0.;
// ==================================================================
#include "truth_reco_var.hh"

int main(int argc, char* argv[]) {
  if (argc==1) {
    cout << "usage: " << argv[0] << " *.root [?i{pf}b]" << endl;
    return 1;
  }
  const std::array<double,2> myy_range{105e3,160e3};

  std::vector<std::unique_ptr<TFile>> mxaods;
  mxaods.reserve(argc-2);
  bool lumi_arg = false;
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
      cout << "\033[36mData\033[0m: " << arg << endl;
      if (a==1+lumi_arg) {
        is_mc = false;
      } else if (is_mc) {
        cerr << "all input files must be either data or MC" << endl;
        return 1;
      }
      const double flumi = std::stod(match[2]);
      lumi += flumi;
      cout << "\033[36mLumi\033[0m: " << flumi << " ipb" << endl;
      mxaods.emplace_back(new TFile(arg,"read"));
      if (mxaods.back()->IsZombie()) return 1;
    } else if (std::regex_search(arg,end,match,mc_re)) { // MC
      cout << "\033[36mMC\033[0m: " << arg << endl;
      if (a==1+lumi_arg) {
        is_mc = true;
      } else if (!is_mc) {
        cerr << "all input files must be either data or MC" << endl;
        return 1;
      }
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
  if (is_mc && !lumi_arg) {
    cerr << "arg error: must specify luminosity for MC" << endl;
    return 1;
  }
  cout << "\n\033[36mLuminosity " << lumi << " ipb\033[0m" << endl << endl;

  // output file for histograms
  auto fout = std::make_unique<TFile>("hists.root","recreate");

  // Histogram definitions ==========================================
#define h_(NAME) TH1D *h_##NAME = new \
  TH1D(#NAME,"",sizeof(b_##NAME)/sizeof(double)-1,b_##NAME);
#define h_truth_(NAME) TH1D *h_##NAME##_truth = new \
  TH1D(#NAME"_truth","",sizeof(b_##NAME)/sizeof(double)-1,b_##NAME);

  double b_HT[] = { 0,30,75,140,200,500 };
  double b_pT_yy[] = { 0,10,15,20,30,45,60,80,100,120,155,200,260,400 };
  double b_pT_j1[] = { 30,40,55,75,95,120,170,400 };
  double b_pT_j2[] = { 30,40,55,95 };
  double b_pT_j3[] = { 30,95 };

  double b_xH[] = { 0,0.25,0.375,0.5,1 };
  double b_x1[] = { 0,0.25,0.375,0.5,1 };
  double b_x2[] = { 0,0.25,0.375,0.5,1 };
  double b_x3[] = { 0,0.25,0.375,0.5,1 };

  h_(HT)
  h_(pT_yy) h_(pT_j1) h_(pT_j2) h_(pT_j3)
  h_(xH) h_(x1) h_(x2) h_(x3)

  // no truth histograms for data
  if (!is_mc) TH1::AddDirectory(false);

  h_truth_(HT)
  h_truth_(pT_yy) h_truth_(pT_j1) h_truth_(pT_j2) h_truth_(pT_j3)
  h_truth_(xH) h_truth_(x1) h_truth_(x2) h_truth_(x3)

  // unsigned nH = 0, n1 = 0, n2 = 0, n3 = 0,
  //   nH_truth = 0, n1_truth = 0, n2_truth = 0, n3_truth = 0;
  unsigned n0 = 0, nn0 = 0;

  for (auto& file : mxaods) { // loop over MxAODs
    cout << "\033[36mMC\033[0m: " << file->GetName() << endl;

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

    // read variables ===============================================
    TTreeReader reader("CollectionTree",file.get());
    TTreeReaderValue<Char_t> _isPassed(reader,"HGamEventInfoAuxDyn.isPassed");
    optional<TTreeReaderValue<Float_t>> _cs_br_fe, _weight;
    optional<TTreeReaderValue<Char_t>> _isFiducial;
    if (is_mc) {
      _cs_br_fe.emplace(reader,"HGamEventInfoAuxDyn.crossSectionBRfilterEff");
      _weight.emplace(reader,"HGamEventInfoAuxDyn.weight");
      _isFiducial.emplace(reader,"HGamTruthEventInfoAuxDyn.isFiducial");
    }

#define VAR_GEN_(NAME, TYPE, STR) \
  var<TTreeReaderValue<TYPE>> _##NAME(reader, STR);
#define VAR_(NAME) VAR_GEN_(NAME, Float_t, #NAME)
#define VAR30_(NAME) VAR_GEN_(NAME, Float_t, #NAME "_30")

    VAR_GEN_(N_j, Int_t, "N_j_30")

    VAR_(m_yy) VAR_(pT_yy)
    VAR30_(pT_j1) VAR30_(pT_j2) VAR30_(pT_j3)
    VAR30_(HT)

    double weight = 1;
    bool is_fiducial = false;

    // loop over events =============================================
    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {

      // selection cut
      if (!*_isPassed) continue;

      // diphoton mass cut
      const auto m_yy = *_m_yy;
      if (!in(m_yy.det,myy_range)) continue;

      if (is_mc) {
        // weight = (**_weight)*(**_cs_br_fe)*n_all_inv*lumi;
        is_fiducial = **_isFiducial && in(m_yy.truth,myy_range);
      }

      const auto nj = *_N_j;

      // FILL HISTOGRAMS ============================================

      const auto HT = _HT/1e3;
      const auto pT_yy = _pT_yy/1e3;
      const auto xH = pT_yy/HT;

      // if (xH > 1) {
      //   test( xH.det )
      //   test( HT.det )
      //   test( pT_yy.det )
      // }
      // if (xH.det > 1) ++nH;
      if (xH.det > 1) {
        if (HT.det == 0) ++n0;
        else ++nn0;
      }

      h_HT->Fill(HT.det,weight);
      h_pT_yy->Fill(pT_yy.det,weight);
      h_xH->Fill(xH.det,weight);

      if (is_fiducial) {
        h_HT_truth->Fill(HT.truth,weight);
        h_pT_yy_truth->Fill(pT_yy.truth,weight);
        h_xH_truth->Fill(xH.truth,weight);
        // if (xH.truth > 1) ++nH_truth;
      }

      if (nj < 1) continue; // 1 jet --------------------------------

      const auto pT_j1 = _pT_j1/1e3;
      const auto x1 = pT_j1/HT;

      // if (x1 > 1) {
      //   test( x1.det )
      //   test( HT.det )
      //   test( pT_j1.det )
      // }
      // if (x1.det > 1) ++n1;

      h_pT_j1->Fill(pT_j1.det,weight);
      h_x1->Fill(x1.det,weight);

      if (is_fiducial && nj.truth >= 1) {
        h_pT_j1_truth->Fill(pT_j1.truth,weight);
        h_x1_truth->Fill(x1.truth,weight);
        // if (x1.truth > 1) ++n1_truth;
      }

      if (nj < 2) continue; // 2 jet --------------------------------

      const auto pT_j2 = _pT_j2/1e3;
      const auto x2 = pT_j2/HT;

      // if (x2.det > 1) ++n2;

      h_pT_j2->Fill(pT_j2.det,weight);
      h_x2->Fill(x2.det,weight);

      if (is_fiducial && nj.truth >= 2) {
        h_pT_j2_truth->Fill(pT_j2.truth,weight);
        h_x2_truth->Fill(x2.truth,weight);
        // if (x2.truth > 1) ++n2_truth;
      }

      if (nj < 3) continue; // 3 jet --------------------------------

      const auto pT_j3 = _pT_j3/1e3;
      const auto x3 = pT_j3/HT;

      // if (x3.det > 1) ++n3;

      h_pT_j3->Fill(pT_j3.det,weight);
      h_x3->Fill(x3.det,weight);

      if (is_fiducial && nj.truth >= 3) {
        h_pT_j3_truth->Fill(pT_j3.truth,weight);
        h_x3_truth->Fill(x3.truth,weight);
        // if (x3.truth > 1) ++n3_truth;
      }

    }
  }

  fout->Write();

  test(n0)
  test(nn0)

  // test( nH )
  // test( nH_truth )
  // test( n1 )
  // test( n1_truth )
  // test( n2 )
  // test( n2_truth )
  // test( n3 )
  // test( n3_truth )
}

