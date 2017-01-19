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

struct file_type_lumi { bool is_mc; double lumi; };
file_type_lumi validate_file_name(const char* fname) {
  static const std::regex data_re(".*/data.*_(\\d*)ipb.*\\.root$");
  static const std::regex mc_re(".*/mc.*\\.root$");
  std::cmatch match;

  const char *begin = fname, *end = fname+std::strlen(fname);
  if (std::regex_search(begin,end,match,data_re)) {
    return { false, std::stod(match[1]) };
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

  re_axes ra("hgam.bins");
// #define a_(name) auto a_##name = ra[#name];
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  hist<ivanp::index_axis<>>
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

#define VAR_STR_(NAME, STR) \
  TTreeReaderValue<Float_t> _##NAME(reader,"HGamEventInfoAuxDyn." STR);
#define VAR_(NAME) VAR_STR_(NAME, #NAME)
#define VAR30_(NAME) VAR_STR_(NAME, #NAME "_30")

    VAR_(m_yy) VAR_(pT_yy) VAR_(yAbs_yy) VAR_(cosTS_yy) VAR_(pTt_yy)
    VAR_(Dy_y_y)

    VAR30_(HT)
    VAR30_(N_j)
    VAR30_(pT_j1)      VAR30_(pT_j2)      VAR30_(pT_j3)
    VAR30_(yAbs_j1)    VAR30_(yAbs_j2)
    VAR30_(Dphi_j_j)   VAR_STR_(Dphi_j_j_signed,"Dphi_j_j_30_signed")
    VAR30_(Dy_j_j)     VAR30_(m_jj)
    VAR30_(sumTau_yyj) VAR30_(maxTau_yyj)
    VAR30_(pT_yyjj)    VAR30_(Dphi_yy_jj)

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

      h_total(0u);

      h_pT_yy   ( pT_yy );
      h_yAbs_yy ( yAbs_yy );
      h_cosTS_yy( cosTS_yy );

      h_Dy_y_y( Dy_y_y );
      h_pTt_yy( *_pTt_yy/1e3 );
      // h_cosTS_pT_yy( cosTS_yy, pT_yy );

      h_N_j_excl( nj );
      h_N_j_incl( 0u );

      h_HT( *_HT/1e3 );

      if (nj == 0) h_pT_yy_0j( pT_yy );

      if (nj < 1) continue; // 1 jet -------------------------------

      h_N_j_incl( 1u );
      
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

      h_N_j_incl( 2u );

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
        if (pT_j3 < 30.) h_VBF(0u);
        if (pT_j3 < 25.) h_VBF(1u);
      }
      if ( (m_jj > 400.) && (dy_jj > 2.8) ) {
        if (pT_j3 < 30.) h_VBF(2u);
      }
      // ------------------------------------------------------------

      if (nj < 3) continue; // 3 jets ------------------------------

      h_N_j_incl( 3u );

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
    for (const auto& h : hist<ivanp::index_axis<>>::all) merge_tmp(h->bins());
    for (const auto& h : re_hist<1>::all) merge_tmp(h->bins());

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
  for (const auto& h : hist<ivanp::index_axis<>>::all) {
    signif(h->bins());
    cout << h << endl;
  }
  for (const auto& h : re_hist<1>::all) {
    signif(h->bins());
    cout << h << endl;
  }

  return 0;
}
