from __future__ import division, print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from root_numpy import root2array, list_trees, array2hist
import pdb
from Utils.Utils import tree_to_df
from ROOT import TH1D
import ROOT
import argparse

parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-i', '--input', nargs='+',
    help='specify input root file'
)
parser.add_argument(
    '-o', '--output', nargs=1, default='./',
    help='specify output dir'
)
args = parser.parse_args()

inputs = args.input
output = args.output[0]

ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/calculateDataEfficiency.C")
ROOT.gROOT.SetBatch(True)

if not os.path.isdir(output):
    os.mkdir(output)

if isinstance(inputs, (list,)):
    treeName = list_trees(inputs[0])
else:
    treeName = list_trees(inputs)
    inputs = [inputs, ]

if (len(treeName) > 1):
    print("more then one tree in file ... specify, which tree to use")
    exit()

# recorded lumi in pb^-1
# from command:$ brilcalc lumi -c web -i /eos/home-d/dwalter/CMSSW_10_6_4/src/ZCounting/TnPPairTreeProducer/production/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU.txt --hltpath="HLT_HIMu17_v*"
Lumi = 211.9

MassBin = 50
MassMin = 66.
MassMax = 116.

ZBBRate = 0.077904
ZBERate = 0.117200
ZEERate = 0.105541

# acceptance selection
selection = None  # 'dilepMass > 66 ' \
# '& dilepMass < 116 ' \

# specify which branches to load
branches = ['nPV', 'dilepMass',  # 'eventNumber', 'run', 'ls',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            'eta1', 'eta2',
            ]

print(">>> Load Events")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
print(">>> Concatenate")
df = pd.concat(df)

hZReco = TH1D("hZReco", "th1d Z reco", MassBin, MassMin, MassMax)

hPassHLTB = TH1D("hPassHLTB", "th1d pass barrel HLT", MassBin, MassMin, MassMax)
hPassHLTE = TH1D("hPassHLTE", "th1d pass endcap HLT", MassBin, MassMin, MassMax)
hFailHLTB = TH1D("hFailHLTB", "th1d fail barrel HLT", MassBin, MassMin, MassMax)
hFailHLTE = TH1D("hFailHLTE", "th1d fail endcap HLT", MassBin, MassMin, MassMax)

hPassSelB = TH1D("hPassSelB", "th1d pass barrel Sel", MassBin, MassMin, MassMax)
hPassSelE = TH1D("hPassSelE", "th1d pass endcap Sel", MassBin, MassMin, MassMax)
hFailSelB = TH1D("hFailSelB", "th1d fail barrel Sel", MassBin, MassMin, MassMax)
hFailSelE = TH1D("hFailSelE", "th1d fail endcap Sel", MassBin, MassMin, MassMax)

hPassGloB = TH1D("hPassGlolB", "th1d pass barrel Glol", MassBin, MassMin, MassMax)
hPassGloE = TH1D("hPassGlolE", "th1d pass endcap Glol", MassBin, MassMin, MassMax)
hFailGloB = TH1D("hFailGlolB", "th1d fail barrel Glol", MassBin, MassMin, MassMax)
hFailGloE = TH1D("hFailGlolE", "th1d fail endcap Glol", MassBin, MassMin, MassMax)

print(">>> Fill Hists")
# --- barrel

for t in df.query("is2HLT==1 & abs(eta1) < 0.9")['dilepMass']:
    hPassHLTB.Fill(t)
    hPassSelB.Fill(t)
    hPassGloB.Fill(t)

for p in df.query("is2HLT==1 & abs(eta2) < 0.9")['dilepMass']:
    hPassHLTB.Fill(p)
    hPassSelB.Fill(p)
    hPassGloB.Fill(p)
    hZReco.Fill(p)

for p in df.query("isSel==1 & abs(eta2) < 0.9")['dilepMass']:
    hFailHLTB.Fill(p)
    hPassSelB.Fill(p)
    hPassGloB.Fill(p)
    hZReco.Fill(p)

for p in df.query("isGlo==1 & abs(eta2) < 0.9")['dilepMass']:
    hFailSelB.Fill(p)
    hPassGloB.Fill(p)

for p in df.query("(isSta==1 | isTrk==1) & abs(eta2) < 0.9")['dilepMass']:
    hFailGloB.Fill(p)

# --- endcap

for t in df.query("is2HLT==1 & abs(eta1) > 0.9")['dilepMass']:
    hPassHLTE.Fill(t)
    hPassSelE.Fill(t)
    hPassGloE.Fill(t)

for p in df.query("is2HLT==1 & abs(eta2) > 0.9")['dilepMass']:
    hPassHLTE.Fill(p)
    hPassSelE.Fill(p)
    hPassGloE.Fill(p)
    hZReco.Fill(p)

for p in df.query("isSel==1 & abs(eta2) > 0.9")['dilepMass']:
    hFailHLTE.Fill(p)
    hPassSelE.Fill(p)
    hPassGloE.Fill(p)
    hZReco.Fill(p)

for p in df.query("isGlo==1 & abs(eta2) > 0.9")['dilepMass']:
    hFailSelE.Fill(p)
    hPassGloE.Fill(p)

for p in df.query("(isSta==1 | isTrk==1) & abs(eta2) > 0.9")['dilepMass']:
    hFailGloE.Fill(p)

print(">>> extract efficiencies")
# --- efficiency extraction

effHLTB = ROOT.calculateDataEfficiency(hPassHLTB, hFailHLTB, output, 0, 4., "HLT", 0, 1, 5, 1, 5, Lumi)
effHLTE = ROOT.calculateDataEfficiency(hPassHLTE, hFailHLTE, output, 0, 4., "HLT", 1, 1, 5, 1, 5, Lumi)

effSelB = ROOT.calculateDataEfficiency(hPassSelB, hFailSelB, output, 0, 4., "Sel", 0, 1, 5, 1, 2, Lumi)
effSelE = ROOT.calculateDataEfficiency(hPassSelE, hFailSelE, output, 0, 4., "Sel", 1, 1, 5, 1, 2, Lumi)

effGloB = ROOT.calculateDataEfficiency(hPassGloB, hFailGloB, output, 0, 4., "Glo", 0, 1, 5, 1, 2, Lumi)
effGloE = ROOT.calculateDataEfficiency(hPassGloE, hFailGloE, output, 0, 4., "Glo", 1, 1, 5, 1, 2, Lumi)

ZBBEff = (effGloB[0] * effGloB[0] * effSelB[0] * effSelB[0] * (1 - (1 - effHLTB[0]) * (1 - effHLTB[0])))
ZBEEff = (effGloB[0] * effGloE[0] * effSelB[0] * effSelE[0] * (1 - (1 - effHLTB[0]) * (1 - effHLTE[0])))
ZEEEff = (effGloE[0] * effGloE[0] * effSelE[0] * effSelE[0] * (1 - (1 - effHLTE[0]) * (1 - effHLTE[0])))

# Statistic Uncertainties (low,high) error propagation
ZBBEff_EStat = [0., 0.]
ZBEEff_EStat = [0., 0.]
ZEEEff_EStat = [0., 0.]
for i in (1, 2):
    ZBBEff_EStat[i - 1] = 2 * ZBBEff * np.sqrt(
        (effGloB[i] / effGloB[0]) ** 2 +
        (effSelB[i] / effSelB[0]) ** 2 +
        ((1 - effHLTB[0]) / (1 - (1 - effHLTE[0]) ** 2) * effHLTB[i]) ** 2
    )
    ZEEEff_EStat[i - 1] = 2 * ZEEEff * np.sqrt(
        (effGloE[i] / effGloE[0]) ** 2 +
        (effSelE[i] / effSelE[0]) ** 2 +
        ((1 - effHLTE[0]) / (1 - (1 - effHLTE[0]) ** 2) * effHLTE[i]) ** 2
    )
    ZBEEff_EStat[i - 1] = ZBEEff * np.sqrt(
        (effGloB[i] / effGloB[0]) ** 2 +
        (effGloE[i] / effGloE[0]) ** 2 +
        (effSelB[i] / effSelB[0]) ** 2 +
        (effSelE[i] / effSelE[0]) ** 2 +
        ((1 - effHLTE[0]) / (1 - (1 - effHLTB[0]) * (1 - effHLTE[0])) * effHLTB[i]) ** 2 +
        ((1 - effHLTB[0]) / (1 - (1 - effHLTB[0]) * (1 - effHLTE[0])) * effHLTE[i]) ** 2
    )

ZEff = (ZBBEff * ZBBRate + ZBEEff * ZBERate + ZEEEff * ZEERate) / (ZBBRate + ZBERate + ZEERate)

ZEff_EStat = [0., 0.]
for i in (0, 1):
    ZEff_EStat[i] = 1./(ZBBRate + ZBERate + ZEERate) * np.sqrt(
        (ZBBRate*ZBBEff_EStat[i])**2 + (ZBERate*ZBEEff_EStat[i])**2 + (ZEERate*ZEEEff_EStat[i])**2
        )

NZReco = hZReco.GetEntries() * 0.99  # assume 1% fake
NZReco_EStat = np.sqrt(hZReco.GetEntries()) * 0.99

NZDeliv = NZReco / ZEff
NZDeliv_EStat = [0., 0.]
for i in (0, 1):
    NZDeliv_EStat[i] = NZDeliv * np.sqrt((NZReco_EStat/NZReco)**2 + (ZEff_EStat[i]/ZEff)**2)

Lumi_Erel = 0.017

ZFid = NZDeliv / Lumi
ZFid_EStat = NZDeliv_EStat[1]/Lumi
ZFid_ELumi = ZFid * Lumi_Erel

print("reconstructed Zs: ", NZReco, " +- ", NZReco_EStat, "(stat.)")
print("Z reconstruction efficiency: ", ZEff, " + ", ZEff_EStat[1], "-", ZEff_EStat[0], " (stat.)")
print("delivered Zs: ", NZDeliv, " + ", NZDeliv_EStat[1], "-", NZDeliv_EStat[0], " (stat.)")
print("Z fid cross section: ", NZDeliv / Lumi, " +- ", ZFid_EStat, "(stat.) +-", ZFid_ELumi, "(Lumi.)")
