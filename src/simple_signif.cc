#include <iostream>
#include <stdexcept>
#include <regex>

#include <boost/optional.hpp>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TKey.h>

using std::cout;
using std::cerr;
using std::endl;
using boost::optional;

template <typename T, typename X>
inline bool in(X x, const std::array<T,2>& a) noexcept {
  return (a[0] < x && x < a[1]);
}

template <typename T>
inline T len(const std::array<T,2>& a) noexcept {
  return (std::get<1>(a) - std::get<0>(a));
}

int main(int argc, const char* argv[]) {
  const std::array<double,2> myy_range{105e3,160e3}, myy_window{121e3,129e3};

  const double data_factor = len(myy_window)/(len(myy_range)-len(myy_window));
  double mc_factor;
  double lumi;

  std::vector<std::pair<std::unique_ptr<TFile>,bool>> mxaods;
  mxaods.reserve(argc-1);

  for (int a=1; a<argc; ++a) { // loop over arguments
    // validate args and parse names of input files
    static const std::regex data_re(
      "^(.*/)?data.*_(\\d*)ipb.*\\.root$", std::regex::optimize);
    static const std::regex mc_re(
      "^(.*/)?mc.*\\.root$", std::regex::optimize);
    std::cmatch match;

    const char *arg = argv[a], *end = arg+std::strlen(arg);
    if (std::regex_search(arg,end,match,data_re)) { // Data
      const double flumi = std::stod(match[2]);
      lumi += flumi;
      cout << "\033[36mData\033[0m: " << arg << endl;
      cout << "\033[36mLumi\033[0m: " << flumi << " ipb" << endl;
      mxaods.emplace_back(std::make_unique<TFile>(arg,"read"),false);
      if (mxaods.back().first->IsZombie())
        throw std::runtime_error("cannot open root file");
    } else if (std::regex_search(arg,end,match,mc_re)) { // MC
      cout << "\033[36mMC\033[0m: " << arg << endl;
      mxaods.emplace_back(std::make_unique<TFile>(arg,"read"),true);
      if (mxaods.back().first->IsZombie())
        throw std::runtime_error("cannot open root file");
    } else {
      cerr << "arg error: unrecognized argument: " << arg << endl;
      return 1;
    }
  }
  if (!mxaods.size()) {
    cerr << "Must specify at least 1 .root file" << endl;
    return 1;
  }
  cout << "\n\033[36mTotal data lumi\033[0m: " << lumi << " ipb" << endl; 

  double sig = 0, bkg = 0;

  for (auto& file : mxaods) { // loop over MxAODs
    const bool is_mc = file.second;
    cout << "\033[36m" << (is_mc ? "MC" : "Data") << "\033[0m: "
         << file.first->GetName() << endl;

    if (is_mc) { // MC
      TIter next(file.first->GetListOfKeys());
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
    }

    // read variables ===============================================
    TTreeReader reader("CollectionTree",file.first.get());
    optional<TTreeReaderValue<Float_t>> _cs_br_fe, _weight;
    if (is_mc) {
      _cs_br_fe.emplace(reader,"HGamEventInfoAuxDyn.crossSectionBRfilterEff");
      _weight.emplace(reader,"HGamEventInfoAuxDyn.weight");
    }
    TTreeReaderValue<Char_t> _isPassed(reader,"HGamEventInfoAuxDyn.isPassed");
    TTreeReaderValue<Float_t> _m_yy(reader,"HGamEventInfoAuxDyn.m_yy");

    while (reader.Next()) { // event loop

      // selection cut
      if (!*_isPassed) continue;

      // diphoton mass cut
      const auto m_yy = *_m_yy;
      if (!in(m_yy,myy_range)) continue;

      const bool is_in_window = in(m_yy,myy_window);

      if (is_mc) { // signal from MC
        if (!is_in_window) continue;
        sig += (**_weight) * (**_cs_br_fe) * mc_factor;
      } else { // background from data
        if (is_in_window) continue;
        bkg += data_factor;
      }

    } // end event loop

  } // end file loop

  cout << endl;
  cout << "signal = " << sig << endl;
  cout << "background = " << bkg << endl;
  cout << "significance = " << sig/std::sqrt(sig+bkg) << endl;
  cout << endl;

} // end main

