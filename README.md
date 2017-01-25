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

# Comments on the code

The [outer loop][L193] in [`main`][L161] loops over the input files. Whether
the file is data or MC is determined from the file name.
If the file is data, the respective luminosity is obtained from the file name
as well.
For MC files, the number of weighted events for the production process is taken
from the `CutFlow_%s_noDalitz_weighted` histogram.

Events in the `CollectionTree` are read using
[`TTreeReader`](https://root.cern.ch/doc/master/classTTreeReader.html)
and
[`TTreeReaderValue`](https://root.cern.ch/doc/master/classTTreeReaderValue.html).

To keep truth and reco level variables together, the [`var`][L121] struct
template is defined. A template specialization is provided for
[`var<TTreeReaderValue<T>>`][L128] for convenient variable definition. It ties
the two value readers for the `"HGamEventInfoAuxDyn."` and
`"HGamTruthEventInfoAuxDyn."` variables together.  The read values are returned
with either `operator*` or `operator()` as a [`var<T>`][L121] struct. The type
`T` is converted to `double` for float point values.  The `operator()` takes a
function as an argument, which is applied to both detector and truth level
values.

[`VAR_`][L243] and [`VAR30_`][L243] macros are defined for convenient
construction of [`var<TTreeReaderValue<T>>`][L128] objects.

The [`fill`][L150] function template is defined to fill histograms with the
values in [`var<T>`][L121].

The truth bin value is filled only if the event passes
`HGamTruthEventInfoAuxDyn.isFiducial`
and if the truth value corresponds to the same bin as at the detector level.

The [`ivanp::binner`](src/binner.hh) class template is used for binning the
variables. The [`hist`][L86], [`re_axis`][L89], and [`re_hist`][L91] aliases
are defined for convenience.
Macro [`h_`][L169] is used for histogram definitions.

struct [`hist_bin`][L33] defines the contents for every bin.
[`operator++`][L47] and [`operator()`][L51] are called when the bin is filled.

[L33]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L33
[L47]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L47
[L51]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L51
[L86]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L86
[L89]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L89
[L91]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L91
[L121]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L121
[L128]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L128
[L150]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L150
[L161]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L161
[L169]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L169
[L193]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L193
[L243]: https://github.com/ivankp/signif3/blob/3adbe10abc6832257898a604bf098f47dfcea098/src/signif.cc#L243
