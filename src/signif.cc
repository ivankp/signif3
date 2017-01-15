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
#include "catstr.hh"
#include "timed_counter.hh"
#include "array_ops.hh"
#include "exception.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

struct file_type_lumi { bool is_mc; double lumi; };
file_type_lumi validate_file_name(const char* fname) {
  static const std::regex data_re(".*/data.*_(\\d*)ipb.*\\.root$");
  static const std::regex mc_re(".*/mc.*\\.root$");
  std::cmatch match;

  const char *begin = fname, *end = fname+std::strlen(fname);
  if (std::regex_search(begin,end,match,data_re)) {
    return { false, atof(match[1].str().c_str()) };
  } else if (std::regex_search(begin,end,match,mc_re)) {
    return { true, 0 };
  } else throw ivanp::exception("Unexpected file name: ",fname);
}

struct hist_bin {
  static double weight;
  double tmp, sig, bkg, signif, purity;
  
  hist_bin(): tmp(0), sig(0), bkg(0), signif(0), purity(0) { }
  inline hist_bin& operator++() noexcept {
    tmp += weight;
    return *this;
  };
};
double hist_bin::weight;

template <typename T>
using axis = ivanp::container_axis<std::vector<T>>;
template <typename... Axes>
using hist_t = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;
template <typename T>
using hist = hist_t<axis<T>>;

template <typename T>
std::ostream& operator<<(std::ostream& o,
  const std::pair<std::string,hist<T>*>& nh
) {
  const auto prec = o.precision();
  const std::ios::fmtflags f( o.flags() );
  o << "\033[32m" << nh.first << "\033[0m\n";
  const auto& h = *nh.second;
  const auto& a = h.axis();
  for (unsigned i=1, n=a.nbins()+2; i<n; ++i) {
    const auto& b = h.bin(i);
    o << "\033[35m[" << a.lower(i) << ',';
    if (i==n-1) o << "∞";
    else o << a.upper(i);
    o << ")\033[0m "
      << round(b.sig) << ' '
      << round(b.bkg) << ' '
      << std::setprecision(2) << std::fixed
      << b.signif << ' '
      << 100*b.purity << '%'
      << std::setprecision(prec) << endl;
    o.flags( f );
  }
  return o;
}

int main(int argc, const char* argv[])
{
  const std::array<double,2>
    mass_range{105e3,160e3},
    mass_window{121e3,129e3};
  const double fw = len(mass_window)/(len(mass_range)-len(mass_window));
  double lumi = 0.;

  hist<int>
    h_total("total",{0,1}),
    h_N_j_excl("N_j_excl",{0,1,2,3,4}),
    h_N_j_incl("N_j_incl",{0,1,2,3,4}),
    h_VBF("VBF",{0,1,2,3});

  hist<double>
    h_pT_yy("pT_yy",{0.,10.,15.,20.,30.,45.,60.,80.,100.,120.,155.,200.,260.,400.}),
    h_yAbs_yy("yAbs_yy",{0.,0.15,0.3,0.45,0.6,0.75,0.9,1.2,1.6,2.4}),
    h_cosTS_yy("cosTS_yy",{0.,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.625,0.75,1.}),
    h_pTt_yy("pTt_yy",{0.,4.,8.,12.,16.,22.,28.,34.,42.,50.,60.,70.,85.,100.,125.,200.}),
    h_Dy_y_y("Dy_y_y",{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,2.}),

    h_HT("HT",{30.,35.,45.,60.,75.,100.,140.,200.,500.}),

    h_pT_j1("pT_j1",{30.,40.,55.,75.,95.,120.,170.,400.}),
    h_pT_j2("pT_j2",{30.,40.,55.,95.}),
    h_pT_j3("pT_j3",{30.,95.}),

    h_yAbs_j1("yAbs_j1",{0.5,1.,1.5,2.,2.5,3.,4.}),
    h_yAbs_j2("yAbs_j2",{0.,1.5,2.5,4.}),

    h_Dphi_j_j("Dphi_j_j",{0., 1.0472, 2.0944, 3.14159}),
    h_Dphi_j_j_signed("Dphi_j_j_signed",{-M_PI,-M_PI_2,0.,M_PI_2,M_PI}),
    h_Dy_j_j("Dy_j_j",{0., 2., 4., 6., 8.8}),
    h_m_jj("m_jj",{0.,200.,5e2,1e3}),

    h_pT_yyjj("pT_yyjj",{0.,20.,40.,150.}),
    h_Dphi_yy_jj("Dphi_yy_jj",{0.,0.05,0.13,0.3,M_PI}),

    h_sumTau_yyj("sumTau_yyj",{0.,2.,4.,8.,12.,17.,25.,40.,80.,150.}),
    h_maxTau_yyj("maxTau_yyj",{0.,1e-3,8.,12.,17.,25.,40.,80.,150.}),

    h_pT_yy_0j("pT_yy_0j",{0.,10.,20.,30.,60.,200.,400.}),
    h_pT_yy_1j("pT_yy_1j",{0.,30.,45.,60.,80.,120.,200.,400.}),
    h_pT_yy_2j("pT_yy_2j",{0.,80.,120.,200.,400.}),
    h_pT_yy_3j("pT_yy_3j",{0.,120.,200.,400.}),

    h_pT_j1_excl("pT_j1_excl",{30.,40.,55.,75.,120.,400.});

    // h_Dphi_Dy_jj("Dphi_Dy_jj",{0.,M_PI_2,M_PI},{0.,2.,8.8}),
    // h_Dphi_pi4_Dy_jj("Dphi_pi4_Dy_jj",{0.,M_PI_2,M_PI},{0.,2.,8.8}),
    // h_cosTS_pT_yy("cosTS_pT_yy",{0.,0.5,1.},{0.,30.,120.,400.}),
    // h_pT_yy_pT_j1("pT_yy_pT_j1",{0.,30.,120.,400.},{30.,65.,400.});

  for (int f=1; f<argc; ++f) {
    const auto file_props = validate_file_name(argv[f]);
    if (file_props.is_mc) {
      cout << "\033[36mMC\033[0m: " << argv[f] << endl;
    } else {
      lumi += file_props.lumi;
      hist_bin::weight = 1;
      cout << "\033[36mData\033[0m: " << argv[f] << endl;
      cout << "\033[36mLumi\033[0m: " << file_props.lumi << " ipb" << endl;
    }

    auto file = std::make_unique<TFile>(argv[f],"read");
    if (file->IsZombie()) return 1;

    double n_all_inv = 1;
    if (file_props.is_mc) {
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
    if (file_props.is_mc) {
      _cs_br_fe.emplace(reader,"HGamEventInfoAuxDyn.crossSectionBRfilterEff");
      _weight .emplace(reader,"HGamEventInfoAuxDyn.weight");
    }
    TTreeReaderValue<Char_t>  _isPassed (reader,"HGamEventInfoAuxDyn.isPassed");

    TTreeReaderValue<Float_t> _m_yy    (reader,"HGamEventInfoAuxDyn.m_yy");
    TTreeReaderValue<Float_t> _pT_yy   (reader,"HGamEventInfoAuxDyn.pT_yy");
    TTreeReaderValue<Float_t> _yAbs_yy (reader,"HGamEventInfoAuxDyn.yAbs_yy");
    TTreeReaderValue<Float_t> _cosTS_yy(reader,"HGamEventInfoAuxDyn.cosTS_yy");
    TTreeReaderValue<Float_t> _pTt_yy  (reader,"HGamEventInfoAuxDyn.pTt_yy");
    TTreeReaderValue<Float_t> _Dy_y_y  (reader,"HGamEventInfoAuxDyn.Dy_y_y");

    TTreeReaderValue<Float_t> _HT(reader,"HGamEventInfoAuxDyn.HT_30");

    TTreeReaderValue<Int_t>   _N_j  (reader,"HGamEventInfoAuxDyn.N_j_30");
    TTreeReaderValue<Float_t> _pT_j1(reader,"HGamEventInfoAuxDyn.pT_j1_30");
    TTreeReaderValue<Float_t> _pT_j2(reader,"HGamEventInfoAuxDyn.pT_j2_30");
    TTreeReaderValue<Float_t> _pT_j3(reader,"HGamEventInfoAuxDyn.pT_j3_30");
    TTreeReaderValue<Float_t> _yAbs_j1(reader,"HGamEventInfoAuxDyn.yAbs_j1_30");
    TTreeReaderValue<Float_t> _yAbs_j2(reader,"HGamEventInfoAuxDyn.yAbs_j2_30");
    TTreeReaderValue<Float_t> _Dphi_j_j(reader,"HGamEventInfoAuxDyn.Dphi_j_j_30");
    TTreeReaderValue<Float_t> _Dphi_j_j_signed(reader,"HGamEventInfoAuxDyn.Dphi_j_j_30_signed");
    TTreeReaderValue<Float_t> _Dy_j_j(reader,"HGamEventInfoAuxDyn.Dy_j_j_30");
    TTreeReaderValue<Float_t> _m_jj  (reader,"HGamEventInfoAuxDyn.m_jj_30");

    TTreeReaderValue<Float_t> _sumTau_yyj(reader,"HGamEventInfoAuxDyn.sumTau_yyj_30");
    TTreeReaderValue<Float_t> _maxTau_yyj(reader,"HGamEventInfoAuxDyn.maxTau_yyj_30");

    TTreeReaderValue<Float_t> _pT_yyjj(reader,"HGamEventInfoAuxDyn.pT_yyjj_30");
    TTreeReaderValue<Float_t> _Dphi_yy_jj(reader,"HGamEventInfoAuxDyn.Dphi_yy_jj_30");

    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
      if (!*_isPassed) continue;
      // is_passed = *_isPassed;

      const double m_yy = *_m_yy;
      if (!in(m_yy,mass_range)) continue;

      if (file_props.is_mc) { // signal from MC
        if (!in(m_yy,mass_window)) continue;
        hist_bin::weight = (**_weight)*(**_cs_br_fe);
      } else { // background from data
        if (in(m_yy,mass_window)) continue;
      }

      // FILL HISTOGRAMS =====================================================

      const unsigned nj = *_N_j;

      const auto pT_yy = *_pT_yy/1e3;
      const auto yAbs_yy = *_yAbs_yy;
      const auto cosTS_yy = std::abs(*_cosTS_yy);
      const auto Dy_y_y = std::abs(*_Dy_y_y);

      h_total(0);

      h_pT_yy   ( pT_yy );
      h_yAbs_yy ( yAbs_yy );
      h_cosTS_yy( cosTS_yy );

      h_Dy_y_y( Dy_y_y );
      h_pTt_yy( *_pTt_yy/1e3 );
      // h_cosTS_pT_yy( cosTS_yy, pT_yy );

      h_N_j_excl( nj );
      h_N_j_incl( 0 );

      h_HT( *_HT/1e3 );

      if (nj == 0) h_pT_yy_0j( pT_yy );

      if (nj < 1) continue; // 1 jet -------------------------------

      h_N_j_incl( 1 );
      
      const auto pT_j1 = *_pT_j1/1e3;
      h_pT_j1( pT_j1 );

      h_yAbs_j1( *_yAbs_j1 );

      h_sumTau_yyj( *_sumTau_yyj/1e3 );
      h_maxTau_yyj( *_maxTau_yyj/1e3 );

      if (nj == 1) {
        h_pT_j1_excl( pT_j1 );
        h_pT_yy_1j( pT_yy );
      }

      // h_pT_yy_pT_j1( pT_yy, pT_j1 );

      if (nj < 2) continue; // 2 jets ------------------------------

      h_N_j_incl( 2 );

      const auto dphi_jj = std::abs(*_Dphi_j_j);
      const auto   dy_jj = std::abs(*_Dy_j_j);
      const auto    m_jj = *_m_jj/1e3;

      h_pT_j2( *_pT_j2/1e3 ); h_yAbs_j2( *_yAbs_j2 );
      
      h_Dphi_yy_jj( M_PI - std::abs(*_Dphi_yy_jj) );

      h_Dphi_j_j_signed( *_Dphi_j_j_signed );
      h_Dphi_j_j( dphi_jj ); // h_HZZ_Dphi_j_j( dphi_jj );
      h_Dy_j_j  (   dy_jj );
      h_m_jj    (    m_jj );

      h_pT_yyjj ( *_pT_yyjj/1e3 );

      if (nj == 2) h_pT_yy_2j( pT_yy );

      // VBF --------------------------------------------------------
      const double pT_j3 = ( nj > 2 ? *_pT_j3/1e3 : 0. );

      if ( (m_jj > 600.) && (dy_jj > 4.0) ) {
        if (pT_j3 < 30.) h_VBF(0);
        if (pT_j3 < 25.) h_VBF(1);
      }
      if ( (m_jj > 400.) && (dy_jj > 2.8) ) {
        if (pT_j3 < 30.) h_VBF(2);
      }
      // ------------------------------------------------------------

      if (nj < 3) continue; // 3 jets ------------------------------

      h_N_j_incl( 3 );

      h_pT_yy_3j( pT_yy );
      h_pT_j3( pT_j3 );

      // h_Dphi_Dy_jj( dphi_jj, dy_jj );
      // h_Dphi_pi4_Dy_jj( phi_pi4(dphi_jj), dy_jj );
    }

    const auto merge_tmp = [=](std::vector<hist_bin>& bins){
      for (auto& b : bins) {
        (file_props.is_mc ? b.sig : b.bkg) += b.tmp*n_all_inv;
        b.tmp = 0;
      }
    };
    for (const auto& h : hist<int>::all) merge_tmp(h.second->bins());
    for (const auto& h : hist<double>::all) merge_tmp(h.second->bins());

  }

  cout << "\n\033[36mTotal lumi\033[0m: " << lumi <<" ipb\n"<< endl;

  const auto signif = [=](std::vector<hist_bin>& bins){
    for (auto& b : bins) {
      b.bkg *= fw;
      b.sig *= lumi;
      const double sb = b.sig + b.bkg;
      b.signif = b.sig/std::sqrt(sb);
      b.purity = b.sig/sb;
    }
  };
  for (const auto& h : hist<int>::all) {
    signif(h.second->bins());
    cout << h << endl;
  }
  for (const auto& h : hist<double>::all) {
    signif(h.second->bins());
    cout << h << endl;
  }

  return 0;
}
