Program for estimating __significance__ and __purity__ of differential
variables for the Hgamma analysis.

# Usage
To compile, just run `make`.
The only dependency is ROOT6.
A C++14 compiler is required. Only tested with gcc 6.2.

The executable is `bin/signif`.
The only arguments to the program are the MxAOD file names.

The variables' binning is specified in the [`hgam.bins`](hgam.bins) file.

# Variables
    m_yy
    pT_yy
    yAbs_yy
    cosTS_yy
    pTt_yy
    Dy_y_y

    N_j_30
    HT_30
    pT_j1_30
    pT_j2_30
    pT_j3_30
    yAbs_j1_30
    yAbs_j2_30
    Dphi_j_j_30
    Dphi_j_j_30_signed
    Dy_j_j_30
    m_jj_30
    sumTau_yyj_30
    maxTau_yyj_30
    pT_yyjj_30
    Dphi_yy_jj_30

# Definitions

**Significance** – _s/√(s+b)_, where _s_ and _b_ are estimates for the numbers
of signal and background events in the signal region (mass window)
respectively.

**Purity** – fraction of signal MC events in a bin (at reco level) for with the
corresponding event at the truth level is in the same bin. The truth level
event must also pass `HGamTruthEventInfoAuxDyn.isFiducial`.

The number of signal events is estimated by counting the number of events with
*m_γγ* within the "mass window" in the Monte Carlo and scaling to the required
luminocity.
The number of MC events is calculated in each bin by
```c++
tmp = Σ ( HGamEventInfoAuxDyn.weight × HGamEventInfoAuxDyn.crossSectionBRfilterEff )
// sum over events for the process

n_all_inv = 1/(total number of weighted events for each production process)
// bin 3 from the CutFlow_%s_noDalitz_weighted histogram

bin = Σ ( tmp * n_all_inv ) // sum over processes
bin *= lumi
```

The number of background events is estimated by counting the number of data
events with *m_γγ* within the "mass range" but not in the "mass window", and
scaling by the factor
```
fw = len(mass_window)/(len(mass_range)-len(mass_window));
```

The "mass range" is 105 to 160 GeV and the "mass window" is 121 to 129 GeV.

# Selection pseudo code
```c++
for ( /* event */ ) {

  if (isPassed) continue;

  if (!(105 < m_yy.det < 160)) continue;

  if (is_mc) {

    if (121 < m_yy.det < 129) {
      // use variable's value for signal estimation for significance
      ++sig;
    }

    if (105 < m_yy.truth < 160) {
      // use variable's value for signal estimation at reco level for purity
      ++reco;
      if ( isFiducial && bin(var.det) == bin(var.truth) ) {
        // use variable's value for signal estimation at truth level for purity
        ++truth;
      }
    }

  } else { // data

    if (!(121 < m_yy.det < 129)) {
      // use variable's value for background estimatio for significancen
      ++bkg;
    }

  }

}

significance = sig/sqrt(sig+bkg);
purity = truth/reco;
```

# Comments on the code

The [outer loop][L308] in [`main`][L239] loops over the input files. Whether
the file is data or MC is [determined][L276] from the file name.
If the file is data, the respective luminosity is obtained from the file name
as well.
For MC files, the number of weighted events for the production process is taken
from the `CutFlow_%s_noDalitz_weighted` histogram.

Events in the `CollectionTree` are read using
[`TTreeReader`](https://root.cern.ch/doc/master/classTTreeReader.html)
and
[`TTreeReaderValue`](https://root.cern.ch/doc/master/classTTreeReaderValue.html).

To keep truth and reco level variables together, the [`var`][L36] struct
template is defined. A template specialization is provided for
[`var<TTreeReaderValue<T>>`][L51] for convenient variable definition. It ties
the two value readers for the `"HGamEventInfoAuxDyn."` and
`"HGamTruthEventInfoAuxDyn."` variables together.  The read values are returned
with either `operator*` or `operator()` as a [`var<T>`][L36] struct. The type
`T` is converted to `double` for float point values.  The `operator()` takes a
function as an argument, which is applied to both detector and truth level
values.

[`VAR_`][L345] and [`VAR30_`][L346] macros are defined for convenient
construction of [`var<TTreeReaderValue<T>>`][L51] objects.

The [`fill`][L196] function templates are defined to fill histograms with the
values in [`var<T>`][L36].

The truth bin value is filled only if the event passes
`HGamTruthEventInfoAuxDyn.isFiducial`
and if the truth value corresponds to the same bin as at the detector level.

The [`ivanp::binner`](src/binner.hh) class template is used for binning the
variables. The [`hist`][L150], [`re_axis`][L153], and [`re_hist`][L155] aliases
are defined for convenience.
Macro [`h_`][L246] is used for histogram definitions.

struct [`hist_bin`][L89] defines the contents for every bin.
The [`operator()`][L102] are called when the bin is filled.

[L89]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L89
[L102]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L51
[L150]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L150
[L153]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L153
[L155]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L155
[L36]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L36
[L51]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L51
[L196]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L196
[L239]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L239
[L276]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L276
[L246]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L246
[L308]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L308
[L345]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L345
[L346]: https://github.com/ivankp/signif3/blob/62cb8deaa15d0c73c2557b4cc014e61d061dff35/src/signif.cc#L346
