from __future__ import division, print_function

import numpy as np
import pandas as pd
import os
from root_numpy import root2array, list_trees, array2hist
import pdb
from Utils.Utils import tree_to_df
from ROOT import TH1D
import ROOT
import json
import glob

import argparse
parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-i', '--input', nargs='+', default=None,
    help='specify input tnp root files'
)
parser.add_argument(
    '--mcCorrections', default='./Resources/MCCorrections.json', type=str,
    help='specify .json file with MC corrections for muon correlations'
)
parser.add_argument(
    '--sigTemplates', default='./Resources/MCTemplates', type=str,
    help='specify directory with MC template root tree files for signal shapes'
)
parser.add_argument(
    '--byLS', default=None, type=str,
    help='specify byLS.json to select specific lumi sections'
)
parser.add_argument(
    '--ptCut', type=float, default=27.,
    help='specify lower pt cut on tag and probe muons'
)
parser.add_argument(
    '--tkIso', default=999., type=float,
    help='specify tracker isolation cut for the HLT criteria'
)
parser.add_argument(
    '--lumi', nargs=3, default=[1, 1, 1], type=float,
    help='specify delivered, recorded luminosity in pb^-1 and its relative uncertainty'
)
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output dir'
)
args = parser.parse_args()

inputs = args.input
if inputs is None:
    inputs = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017H_LowPU_V02/SingleMuon/TnPPairTrees_2017H_LowPU_V02/191014_135140/0000/output_0_*")

#inputs of same sign selection
inputs_ss = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017H_LowPUSameCharge_V02/SingleMuon/TnPPairTrees_2017H_LowPUSameCharge_V02/191020_141054/0000/output_0_*")

#high PU input
inputs_os_highPu2017B = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017B_V02/SingleMuon/TnPPairTrees_2017B_V02/191022_123558/0000/output_0_*")
inputs_os_highPu2017C = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017C_V02/SingleMuon/TnPPairTrees_2017C_V02/191022_124705/0000/output_0_*")
inputs_os_highPu2017D = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017D_V02/SingleMuon/TnPPairTrees_2017D_V02/191022_124751/0000/output_0_*")
inputs_os_highPu2017E = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017E_V02/SingleMuon/TnPPairTrees_2017E_V02/191022_124830/0000/output_0_*")
inputs_os_highPu2017F = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017F_V02/SingleMuon/TnPPairTrees_2017F_V02/191022_124919/0000/output_0_*")

inputs_ss_highPu2017B = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017B_SameCharge_V02_1/SingleMuon/TnPPairTrees_2017B_SameCharge_V02_1/191109_164559/0000/output_0_*")
inputs_ss_highPu2017C = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017C_SameCharge_V02/SingleMuon/TnPPairTrees_2017C_SameCharge_V02/191102_165018/0000/output_0_*")
inputs_ss_highPu2017D = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017D_SameCharge_V02/SingleMuon/TnPPairTrees_2017D_SameCharge_V02/191109_164614/0000/output_0_*")
inputs_ss_highPu2017E = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017E_SameCharge_V02/SingleMuon/TnPPairTrees_2017E_SameCharge_V02/191109_164627/0000/output_0_*")
inputs_ss_highPu2017F = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/TnPPairTrees_2017F_SameCharge_V02/SingleMuon/TnPPairTrees_2017F_SameCharge_V02/191109_164641/0000/output_0_*")

inputs_ss_highPu = [inputs_ss_highPu2017B, inputs_ss_highPu2017C, inputs_ss_highPu2017D, inputs_ss_highPu2017E, inputs_ss_highPu2017F]
inputs_os_highPu = [inputs_os_highPu2017B, inputs_os_highPu2017C, inputs_os_highPu2017D, inputs_os_highPu2017E, inputs_os_highPu2017F]

output = args.output

mcCorr = args.mcCorrections
sigTemplates = args.sigTemplates
byLS_file = args.byLS
ptCut = args.ptCut
tkIsoCut = args.tkIso
lumiDel, lumiRec, lumi_Erel = args.lumi

lumiDel, lumiRec, lumi_Erel = 200, 199, 0.017
ptCut = 30
byLS_file = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU.txt"
#sigTemplates_DY = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/MCTemplates_DY50TuneCP1andCP5FlatPU0to75_Pt30HLT_L1SingleMu25/"
sigTemplates_DY = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/MCTemplates_DY50TuneCP1andCP5FlatPU0to75_Pt30HLT_IsoMu27/"

mcCorr = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/MCTemplates_DY50TuneCP1andCP5FlatPU0to75/MCCorrections_nPU_TuneCP5FlatPU0to75_TightIDIsoMu27Pt30/MCCorrections_nPU_ZTightID_IsoMu27.json"


bkgTemplateQCD = ""
bkgTemplateTT = ""



#ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/calculateDataEfficiency.C")
#ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/calculateFakeRate.C")
ROOT.gROOT.SetBatch(True)


if os.path.isdir(output):
    print("output dir already exists, please remove or specify another name")
    exit()
os.mkdir(output)

if isinstance(inputs, (list,)):
    treeName = list_trees(inputs[0])
else:
    treeName = list_trees(inputs)
    inputs = [inputs, ]

if (len(treeName) > 1):
    print("more then one tree in file ... specify, which tree to use")
    exit()

MassBin = 50
MassMin = 66.
MassMax = 116.

ZBBRate = 0.077904
ZBERate = 0.117200
ZEERate = 0.105541

# acceptance selection
selection = 'dilepMass > 66 & dilepMass < 116 & pt1 > {0} & pt2 > {0} & tkIso1/pt1 < {1}'.format(ptCut, tkIsoCut)  # ' ' \
# '  ' \

# specify which branches to load
branches = ['nPV', 'dilepMass',  # 'eventNumber', 'run', 'ls',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            'pt2', 'tkIso2',
            'eta1', 'eta2',
            'run','ls'
            ]

print(">>> Load opposite events in low PU")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
print(">>> Concatenate")
df = pd.concat(df)

print(">>> Load same sign events in low PU")
df_ss = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs_ss]
print(">>> Concatenate")
df_ss = pd.concat(df_ss)

print(">>> Load same sign events in high PU")
list_ss_highPU = []
for input_ss_highPu in inputs_ss_highPu:
    df_ss_highPu = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in input_ss_highPu]
    print(">>> Concatenate")
    list_ss_highPU.append(pd.concat(df_ss_highPu))

print(">>> Load opposite sign events in high PU")
list_os_highPU = []
for input_ss_highPu in inputs_os_highPu:
    df_os_highPu = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in input_ss_highPu]
    print(">>> Concatenate")
    list_os_highPU.append(pd.concat(df_os_highPu))

if byLS_file is not None:
    with open(byLS_file) as json_file:
        byLS = json.load(json_file)
    print(">>> select good ls from byLS.json")
    df = df.query('|'.join(['(run == {0} & ls >= {1} & ls <= {2})'.format(run,ls_start,ls_end) for run, ls_bounds in byLS.items() for ls_start, ls_end in ls_bounds]))
    df_ss = df_ss.query('|'.join(['(run == {0} & ls >= {1} & ls <= {2})'.format(run,ls_start,ls_end) for run, ls_bounds in byLS.items() for ls_start, ls_end in ls_bounds]))


# df['is2HLT'] = df['is2HLT'] * (df['tkIso2'] / df['pt2'] < tkIsoCut)
# df['isSel'] = df['isSel'] + df['is2HLT'] * (df['tkIso2'] / df['pt2'] > tkIsoCut)
#
# df_ss['is2HLT'] = df_ss['is2HLT'] * (df_ss['tkIso2'] / df_ss['pt2'] < tkIsoCut)
# df_ss['isSel'] = df_ss['isSel'] + df_ss['is2HLT'] * (df_ss['tkIso2'] / df_ss['pt2'] > tkIsoCut)
#
# df_ss_highPu['is2HLT'] = df_ss_highPu['is2HLT'] * (df_ss_highPu['tkIso2'] / df_ss_highPu['pt2'] < tkIsoCut)
# df_ss_highPu['isSel'] = df_ss_highPu['isSel'] + df_ss_highPu['is2HLT'] * (df_ss_highPu['tkIso2'] / df_ss_highPu['pt2'] > tkIsoCut)
#
# df_os_highPu['is2HLT'] = df_os_highPu['is2HLT'] * (df_os_highPu['tkIso2'] / df_os_highPu['pt2'] < tkIsoCut)
# df_os_highPu['isSel'] = df_os_highPu['isSel'] + df_os_highPu['is2HLT'] * (df_os_highPu['tkIso2'] / df_os_highPu['pt2'] > tkIsoCut)

hZReco = TH1D("hZReco", "th1d Z reco", MassBin, MassMin, MassMax)
hZReco_ss = TH1D("hZReco_ss", "th1d Z reco same sign", MassBin, MassMin, MassMax)
hZReco_ss_highPu = TH1D("hZReco_ss_highPu", "th1d Z reco same sign high pileup", MassBin, MassMin, MassMax)
hZReco_os_highPu = TH1D("hZReco_os_highPu", "th1d Z reco opposite sign high pileup", MassBin, MassMin, MassMax)

# hPV_os_lowPu = TH1D("hPV_data", "th1d number of primary vertices shape", 75, -0.5, 74.5)
# hPV_ss_lowPu = TH1D("hPV_data", "th1d number of primary vertices shape", 75, -0.5, 74.5)
# hPV_os_highPu = TH1D("hPV_data", "th1d number of primary vertices shape", 75, -0.5, 74.5)
# hPV_ss_highPu = TH1D("hPV_data", "th1d number of primary vertices shape", 75, -0.5, 74.5)
#
# hPassHLTB = TH1D("hPassHLTB", "th1d pass barrel HLT", MassBin, MassMin, MassMax)
# hPassHLTE = TH1D("hPassHLTE", "th1d pass endcap HLT", MassBin, MassMin, MassMax)
# hFailHLTB = TH1D("hFailHLTB", "th1d fail barrel HLT", MassBin, MassMin, MassMax)
# hFailHLTE = TH1D("hFailHLTE", "th1d fail endcap HLT", MassBin, MassMin, MassMax)
#
# hPassSelB = TH1D("hPassSelB", "th1d pass barrel Sel", MassBin, MassMin, MassMax)
# hPassSelE = TH1D("hPassSelE", "th1d pass endcap Sel", MassBin, MassMin, MassMax)
# hFailSelB = TH1D("hFailSelB", "th1d fail barrel Sel", MassBin, MassMin, MassMax)
# hFailSelE = TH1D("hFailSelE", "th1d fail endcap Sel", MassBin, MassMin, MassMax)
#
# hPassGloB = TH1D("hPassGlolB", "th1d pass barrel Glol", MassBin, MassMin, MassMax)
# hPassGloE = TH1D("hPassGlolE", "th1d pass endcap Glol", MassBin, MassMin, MassMax)
# hFailGloB = TH1D("hFailGlolB", "th1d fail barrel Glol", MassBin, MassMin, MassMax)
# hFailGloE = TH1D("hFailGlolE", "th1d fail endcap Glol", MassBin, MassMin, MassMax)
#
# print(">>> Fill Hists")
# for nPV in df['nPV']:
#     hPV_os_lowPu.Fill(nPV)
#
# hPV_os_lowPu.Scale(1. / hPV_os_lowPu.Integral())

# fill tag and probe histograms
# --- barrel
# for t in df.query("is2HLT==1 & abs(eta1) < 0.9")['dilepMass']:
#     hPassHLTB.Fill(t)
#     hPassSelB.Fill(t)
#     hPassGloB.Fill(t)
#
# for p in df.query("is2HLT==1 & abs(eta2) < 0.9")['dilepMass']:
#     hPassHLTB.Fill(p)
#     hPassSelB.Fill(p)
#     hPassGloB.Fill(p)
#
# for p in df.query("isSel==1 & abs(eta2) < 0.9")['dilepMass']:
#     hFailHLTB.Fill(p)
#     hPassSelB.Fill(p)
#     hPassGloB.Fill(p)
#
# for p in df.query("isGlo==1 & abs(eta2) < 0.9")['dilepMass']:
#     hFailSelB.Fill(p)
#     hPassGloB.Fill(p)
#
# for p in df.query("(isSta==1 | isTrk==1) & abs(eta2) < 0.9")['dilepMass']:
#     hFailGloB.Fill(p)
#
# # --- endcap
# for t in df.query("is2HLT==1 & abs(eta1) > 0.9")['dilepMass']:
#     hPassHLTE.Fill(t)
#     hPassSelE.Fill(t)
#     hPassGloE.Fill(t)
#
# for p in df.query("is2HLT==1 & abs(eta2) > 0.9")['dilepMass']:
#     hPassHLTE.Fill(p)
#     hPassSelE.Fill(p)
#     hPassGloE.Fill(p)
#
# for p in df.query("isSel==1 & abs(eta2) > 0.9")['dilepMass']:
#     hFailHLTE.Fill(p)
#     hPassSelE.Fill(p)
#     hPassGloE.Fill(p)
#
# for p in df.query("isGlo==1 & abs(eta2) > 0.9")['dilepMass']:
#     hFailSelE.Fill(p)
#     hPassGloE.Fill(p)
#
# for p in df.query("(isSta==1 | isTrk==1) & abs(eta2) > 0.9")['dilepMass']:
#     hFailGloE.Fill(p)

df_ss_highPu = pd.concat(list_ss_highPU)

for p in df_ss_highPu.query("is2HLT==1 | isSel==1")['dilepMass']:
    hZReco_ss_highPu.Fill(p)


yieldfitter = ROOT.RooFitter(ptCut, ptCut, 2, 5, output, sigTemplates_DY+ "template_ZYield.root", 0)
yieldfitter.fit_backgroundmodel(hZReco_ss_highPu)

# fill z yield same histograms
for p in df.query("is2HLT==1 | isSel==1")['dilepMass']:
    hZReco.Fill(p)

for p in df_ss.query("is2HLT==1 | isSel==1")['dilepMass']:
    hZReco_ss.Fill(p)

res = yieldfitter.fit_simultanious(hZReco, hZReco_ss, "H")
results = {}
results["H"] = pd.DataFrame([[res[i] for i in range(len(res))] + [df.query("is2HLT==1 | isSel==1")['nPV'].mean()]],
    columns=["fr","fr_err","tf","tf_err","chi2_OS","Chi2_SS","avgPU"])

eras = ["B", "C", "D", "E", "F"]
iStep = 100000
# loop over eras B to F
for df_os_highPu, df_ss_highPu, era in zip(list_os_highPU, list_ss_highPU, eras):

    dfReco_os_highPu = df_os_highPu.query("is2HLT==1 | isSel==1").sort_values(by=['nPV'])
    dfReco_ss_highPu = df_ss_highPu.query("is2HLT==1 | isSel==1").sort_values(by=['nPV'])

    result = []
    iMeasure = 0
    iLast = 0
    # loop over measurements
    while (iMeasure+1)*iStep < len(dfReco_os_highPu):
        print("era {0} do measurement {1}/{2}".format(era, iMeasure, len(dfReco_os_highPu)//iStep) )
        print("{0} < {1}".format((iMeasure+1)*iStep, len(dfReco_os_highPu)))



        hZReco_os_highPu = TH1D("hZReco_os_highPu", "th1d Z reco opposite sign high pileup", MassBin, MassMin, MassMax)
        hZReco_ss_highPu = TH1D("hZReco_ss_highPu", "th1d Z reco same sign high pileup", MassBin, MassMin, MassMax)

        selected = dfReco_os_highPu[iMeasure*iStep:(iMeasure+1)*iStep]
        npv_min = selected['nPV'].min()
        npv_max = selected['nPV'].max()
        npv_mean = selected['nPV'].mean()

        npv_cummmean = dfReco_os_highPu[:(iMeasure+1)*iStep]['nPV'].mean()

        selected_ss = dfReco_ss_highPu[dfReco_ss_highPu['nPV'].expanding().mean() < npv_cummmean]['dilepMass']

        for p in selected['dilepMass']:
            hZReco_os_highPu.Fill(p)

        for p in selected_ss[iLast:]:
            hZReco_ss_highPu.Fill(p)

        iLast = len(selected_ss)

        res = yieldfitter.fit_simultanious(hZReco_os_highPu, hZReco_ss_highPu, "{0}_{1}".format(era,iMeasure))

        result.append(([res[i] for i in range(len(res))]+ [npv_mean]))

        iMeasure += 1
        hZReco_os_highPu.Delete()
        hZReco_ss_highPu.Delete()

    results[era] = pd.DataFrame(result, columns=["fr","fr_err","tf","tf_err","chi2_OS","Chi2_SS","nPV"])


import pickle
pickle.dump(results, file(output+"/results.p","w"))

exit()


# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 2, 1, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 3, 1, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 4, 1, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 5, 1, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
#
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 1, 2, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 2, 2, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 3, 2, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 4, 2, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 5, 2, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
#
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 1, 3, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 2, 3, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 3, 3, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 4, 3, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 5, 3, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
#
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 1, 4, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 2, 4, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 3, 4, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 4, 4, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)
# xxx = ROOT.calculateZYield(hZReco, hZReco_ss, hZReco_ss_highPu, 2, 5, 4, ptCut, ptCut, 0, sigTemplates_DY + "template_ZYield.root", output)


print(">>> extract efficiencies")
# --- efficiency extraction

# --- MC corrections
with open(mcCorr) as json_file:
    corr = json.load(json_file)

result = []


def calculateEfficiencies(sigShapes, bkgShapes, outp):
    if not os.path.isdir(outp):
        os.mkdir(outp)

    effGloB = ROOT.calculateDataEfficiency(hPassGloB, hFailGloB, outp, 0, "Glo", 0,
                                           sigShapes[0], bkgShapes[0], sigShapes[1], bkgShapes[1], ptCut, ptCut, 0,
                                           # hPV,
                                           lumiRec, sigTemplates_DY + "template_Glo.root", bkgTemplateQCD, bkgTemplateTT)
    effGloE = ROOT.calculateDataEfficiency(hPassGloE, hFailGloE, outp, 0, "Glo", 1,
                                           sigShapes[0], bkgShapes[0], sigShapes[1], bkgShapes[1], ptCut, ptCut, 0,
                                           # hPV,
                                           lumiRec, sigTemplates_DY + "template_Glo.root", bkgTemplateQCD, bkgTemplateTT)
    effSelB = ROOT.calculateDataEfficiency(hPassSelB, hFailSelB, outp, 0, "Sel", 0,
                                           sigShapes[2], bkgShapes[2], sigShapes[3], bkgShapes[3], ptCut, ptCut, 0,
                                           # hPV,
                                           lumiRec, sigTemplates_DY + "template_Sel.root", bkgTemplateQCD, bkgTemplateTT)
    effSelE = ROOT.calculateDataEfficiency(hPassSelE, hFailSelE, outp, 0, "Sel", 1,
                                           sigShapes[2], bkgShapes[2], sigShapes[3], bkgShapes[3], ptCut, ptCut, 0,
                                           # hPV,
                                           lumiRec, sigTemplates_DY + "template_Sel.root", bkgTemplateQCD, bkgTemplateTT)
    effHLTB = ROOT.calculateDataEfficiency(hPassHLTB, hFailHLTB, outp, 0, "HLT", 0,
                                           sigShapes[4], bkgShapes[4], sigShapes[5], bkgShapes[5], ptCut, ptCut, 0,
                                           # hPV,
                                           lumiRec, sigTemplates_DY + "template_HLT.root", bkgTemplateQCD, bkgTemplateTT)
    effHLTE = ROOT.calculateDataEfficiency(hPassHLTE, hFailHLTE, outp, 0, "HLT", 1,
                                           sigShapes[4], bkgShapes[4], sigShapes[5], bkgShapes[5], ptCut, ptCut, 0,
                                           # hPV,
                                           lumiRec, sigTemplates_DY + "template_HLT.root", bkgTemplateQCD, bkgTemplateTT)

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

    avgPU = df['nPV'].mean()
    ZBBEffMC = ZBBEff - (corr['BB_a'] * avgPU + corr['BB_b'])
    ZBEEffMC = ZBEEff - (corr['BE_a'] * avgPU + corr['BE_b'])
    ZEEEffMC = ZEEEff - (corr['EE_a'] * avgPU + corr['EE_b'])

    ZEff = (ZBBEff * ZBBRate + ZBEEff * ZBERate + ZEEEff * ZEERate) / (ZBBRate + ZBERate + ZEERate)
    ZEffMC = (ZBBEffMC * ZBBRate + ZBEEffMC * ZBERate + ZEEEffMC * ZEERate) / (ZBBRate + ZBERate + ZEERate)

    ZEff_EStat = [0., 0.]
    for i in (0, 1):
        ZEff_EStat[i] = 1. / (ZBBRate + ZBERate + ZEERate) * np.sqrt(
            (ZBBRate * ZBBEff_EStat[i]) ** 2 + (ZBERate * ZBEEff_EStat[i]) ** 2 + (ZEERate * ZEEEff_EStat[i]) ** 2
        )

    NZReco = hZReco.GetEntries() * 0.99  # assume 1% fake
    NZReco_EStat = np.sqrt(hZReco.GetEntries()) * 0.99

    deadtime = lumiRec / lumiDel

    NZDeliv = NZReco / (ZEff * deadtime)
    NZDelivMC = NZReco / (ZEffMC * deadtime)
    NZDelivMC_EStat = [0., 0.]
    for i in (0, 1):
        NZDelivMC_EStat[i] = NZDeliv * np.sqrt((NZReco_EStat / NZReco) ** 2 + (ZEff_EStat[i] / ZEffMC) ** 2)

    ZFid = NZDeliv / lumiDel
    ZFidMC = NZDelivMC / lumiDel
    ZFidMC_EStat = NZDelivMC_EStat[1] / lumiDel
    ZFidMC_ELumi = ZFidMC * lumi_Erel

    result.append([sigShapes, bkgShapes,
                   effGloB[0], effGloB[1], effGloE[0], effGloE[1],
                   effSelB[0], effSelB[1], effSelE[0], effSelE[1],
                   effHLTB[0], effHLTB[1], effHLTE[0], effHLTE[1],
                   effGloB[3], effGloB[4], effGloE[3], effGloE[4],
                   effSelB[3], effSelB[4], effSelE[3], effSelE[4],
                   effHLTB[3], effHLTB[4], effHLTE[3], effHLTE[4],
                   ZBBEff, ZBBEff_EStat, ZBEEff, ZBEEff_EStat, ZEEEff, ZEEEff_EStat,
                   ZEff, ZEff_EStat, ZEffMC,
                   NZReco, NZReco_EStat,
                   NZDeliv, NZDelivMC, NZDelivMC_EStat[1],
                   ZFid, ZFidMC, ZFidMC_EStat, ZFidMC_ELumi
                   ])


#calculateEfficiencies([2, 2, 2, 2, 2, 2], [5, 2, 5, 5, 5, 5], output + "/effNominal")
calculateEfficiencies([2, 2, 2, 2, 2, 2], [2, 2, 1, 1, 5, 5], output + "/effNominal")
#calculateEfficiencies([1, 1, 1, 1, 1, 1], [7, 7, 7, 7, 7, 7], output + "/effVarSig")
#calculateEfficiencies([1, 1, 1, 1, 1, 1], [5, 2, 5, 5, 5, 5], output + "/effVarSig")

results = pd.DataFrame(result,
                       columns=['sigShapes', 'bkgShapes',
                                'effGloB', 'effGloB_E', 'effGloE', 'effGloE_E',
                                'effSelB', 'effSelB_E', 'effSelE', 'effSelE_E',
                                'effHLTB', 'effHLTB_E', 'effHLTE', 'effHLTE_E',
                                'effGloB_chi2Pass', 'effGloB_chi2Fail', 'effGloE_chi2Pass', 'effGloE_chi2Fail',
                                'effSelB_chi2Pass', 'effSelB_chi2Fail', 'effSelE_chi2Pass', 'effSelE_chi2Fail',
                                'effHLTB_chi2Pass', 'effHLTB_chi2Fail', 'effHLTE_chi2Pass', 'effHLTE_chi2Fail',
                                'ZBBEff', 'ZBBEff_EStat', 'ZBEEff', 'ZBEEff_EStat', 'ZEEEff', 'ZEEEff_EStat',
                                'ZEff', 'ZEff_EStat', 'ZEffMC',
                                'NZReco', 'NZReco_EStat',
                                'NZDeliv', 'NZDelivMC', 'NZDelivMC_EStat',
                                'ZFid', 'ZFidMC', 'ZFidMC_EStat', 'ZFidMC_ELumi'
                                ])

with open(output + '/result.csv', 'w') as file:
    results.to_csv(file)
