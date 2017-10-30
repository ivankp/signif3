#!/bin/bash

MxAOD=/eos/atlas/atlascerngroupdisk/phys-higgs/HSG1/MxAOD/h015d

./bin/signif ${1-signif.bins} \
  $(find ${MxAOD}/data{15,16} -maxdepth 1 -type f -name '*.root') \
  ${MxAOD}/mc15c/mc15c.PowhegPy8_NNLOPS_ggH125.MxAODDetailed.p3015.h015d.root \
  ${MxAOD}/mc15c/mc15c.PowhegPy8_NNPDF30_VBFH125.MxAODDetailed.p3015.h015d.root \
  ${MxAOD}/mc15c/mc15c.PowhegPy8_WmH125J.MxAODDetailed.p3015.h015d.root \
  ${MxAOD}/mc15c/mc15c.PowhegPy8_WpH125J.MxAODDetailed.p3015.h015d.root \
  ${MxAOD}/mc15c/mc15c.PowhegPy8_ZH125J.MxAODDetailed.p3015.h015d.root \
  ${MxAOD}/mc15c/mc15c.PowhegPy8_ggZH125.MxAODDetailed.p3015.h015d.root \
  ${MxAOD}/mc15c/mc15c.aMCnloPy8_ttH125.MxAODDetailed.p2908.h015d.root \
  ${MxAOD}/mc15c/mc15c.aMCnloPy8_bbH125_yb2.MxAODDetailed.p2908.h015d.root \
  ${MxAOD}/mc15c/mc15c.aMCnloPy8_bbH125_ybyt.MxAODDetailed.p2908.h015d.root

