#include <iostream>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
  TFile file(argv[1],"read");
  
  TTreeReader reader("CollectionTree",&file);
  TTreeReaderValue<Char_t> _isPassed(reader,"HGamEventInfoAuxDyn.isPassed");
  TTreeReaderValue<Int_t> _N_j(reader,"HGamEventInfoAuxDyn.N_j_30");
  TTreeReaderValue<Float_t> _HT(reader,"HGamEventInfoAuxDyn.HT_30");
  TTreeReaderValue<Float_t> _pT_j1(reader,"HGamEventInfoAuxDyn.pT_j1_30");
  
  cout << std::boolalpha;

  for (auto event = 0ull; reader.Next(); ++event) {
    if (*_N_j < 1) continue;

    const double HT = *_HT*1e-3;
    const double pT_j1 = *_pT_j1*1e-3;

    if (HT < pT_j1) {
      cout << event << ':'
           << " HT=" << HT
           << " pT_j1=" << pT_j1
           << " passed=" << *_isPassed
           << endl;
    }
  }

  return 0;
}
