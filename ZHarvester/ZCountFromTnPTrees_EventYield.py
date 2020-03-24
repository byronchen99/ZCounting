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
ROOT.gROOT.SetBatch(True)

import argparse
parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output dir'
)
args = parser.parse_args()

output=args.output

storage="/pnfs/desy.de/cms/tier2/store/user/dwalter/"

inputs = glob.glob(storage+"TnPPairTrees_2017H_LowPU_V02/SingleMuon/TnPPairTrees_2017H_LowPU_V02/191014_135140/0000/output_0_*")

#inputs of same sign selection
inputs_ss = glob.glob(storage+"TnPPairTrees_2017H_LowPUSameCharge_V02/SingleMuon/TnPPairTrees_2017H_LowPUSameCharge_V02/191020_141054/0000/output_0_*")

#high PU input
inputs_os_highPu2017B = glob.glob(storage+"TnPPairTrees_2017B_V02/SingleMuon/TnPPairTrees_2017B_V02/191022_123558/0000/output_0_*")
inputs_os_highPu2017C = glob.glob(storage+"TnPPairTrees_2017C_V02/SingleMuon/TnPPairTrees_2017C_V02/191022_124705/0000/output_0_*")
inputs_os_highPu2017D = glob.glob(storage+"TnPPairTrees_2017D_V02/SingleMuon/TnPPairTrees_2017D_V02/191022_124751/0000/output_0_*")
inputs_os_highPu2017E = glob.glob(storage+"TnPPairTrees_2017E_V02/SingleMuon/TnPPairTrees_2017E_V02/191022_124830/0000/output_0_*")
inputs_os_highPu2017F = glob.glob(storage+"TnPPairTrees_2017F_V02/SingleMuon/TnPPairTrees_2017F_V02/191022_124919/0000/output_0_*")

inputs_ss_highPu2017B = glob.glob(storage+"TnPPairTrees_2017B_SameCharge_V02_1/SingleMuon/TnPPairTrees_2017B_SameCharge_V02_1/191109_164559/0000/output_0_*")
inputs_ss_highPu2017C = glob.glob(storage+"TnPPairTrees_2017C_SameCharge_V02/SingleMuon/TnPPairTrees_2017C_SameCharge_V02/191102_165018/0000/output_0_*")
inputs_ss_highPu2017D = glob.glob(storage+"TnPPairTrees_2017D_SameCharge_V02/SingleMuon/TnPPairTrees_2017D_SameCharge_V02/191109_164614/0000/output_0_*")
inputs_ss_highPu2017E = glob.glob(storage+"TnPPairTrees_2017E_SameCharge_V02/SingleMuon/TnPPairTrees_2017E_SameCharge_V02/191109_164627/0000/output_0_*")
inputs_ss_highPu2017F = glob.glob(storage+"TnPPairTrees_2017F_SameCharge_V02/SingleMuon/TnPPairTrees_2017F_SameCharge_V02/191109_164641/0000/output_0_*")

inputs_ss_highPu = [inputs_ss_highPu2017B, inputs_ss_highPu2017C, inputs_ss_highPu2017D, inputs_ss_highPu2017E, inputs_ss_highPu2017F]
inputs_os_highPu = [inputs_os_highPu2017B, inputs_os_highPu2017C, inputs_os_highPu2017D, inputs_os_highPu2017E, inputs_os_highPu2017F]

ptCut = 30
byLS_file = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU_lowPU.txt"

sigTemplates_DY_flatPU = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/MCTemplates/DY50_TuneCP5MadgraphMLMPythia8_Pt30_L1SingleMu25/"
sigTemplates_DY_L1SingleMu25 = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/MCTemplates/DY50_TuneCP5amcatnloFXFXPythia8_Pt30_L1SingleMu25/"
sigTemplates_DY = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/MCTemplates/DY50_TuneCP5amcatnloFXFXPythia8_IsoMu27/"

sigTemplates_TT = "/nfs/dust/cms/user/dwalter/data/Lumi/ZMonitoring2017/MCTemplates/TTTo2L2Nu_TuneCP5PowhegPythia8_IsoMu27/"

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

# specify which branches to load
branches = ['nPV', 'dilepMass',  # 'eventNumber', 'run', 'ls',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            'pt2', 'tkIso2',
            'eta1', 'eta2',
            'run','ls'
            ]

selection = 'dilepMass > 66 & dilepMass < 116 & pt1 > {0} & pt2 > {0}'.format(ptCut)

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


hPV = TH1D("hPV_data", "th1d number of primary vertices shape", 75, -0.5, 74.5)

hZReco = TH1D("hZReco", "th1d Z reco", MassBin, MassMin, MassMax)
hZReco_ss = TH1D("hZReco_ss", "th1d Z reco same sign", MassBin, MassMin, MassMax)

df_ss_highPu = pd.concat(list_ss_highPU)
for p in df_ss_highPu.query("is2HLT==1 | isSel==1")['dilepMass']:
    hZReco_ss.Fill(p)

for p in df_ss.query("is2HLT==1 | isSel==1")['dilepMass']:
    hZReco_ss.Fill(p)

yieldfitter = ROOT.RooFitter(ptCut, ptCut, 2, 5, output,
    sigTemplates_DY + "template_ZYield.root")
yieldfitter.fit_backgroundmodel(hZReco_ss)

print(">>> Fill Hists 2017 H")
for p in df.query("is2HLT==1 | isSel==1")['dilepMass']:
    hZReco.Fill(p)
for nPV in df['nPV']:
    hPV.Fill(nPV)

hZReco_ss.Delete()
hZReco_ss = TH1D("hZReco_ss", "th1d Z reco same sign", MassBin, MassMin, MassMax)
for p in df_ss.query("is2HLT==1 | isSel==1")['dilepMass']:
    hZReco_ss.Fill(p)

npv_mean = df.query("is2HLT==1 | isSel==1")['nPV'].mean()
result = []

print("fit flat PU template w/o nPV reweighting")
res = yieldfitter.fit_simultanious(hZReco, hZReco_ss, "H_0")
result.append(([res[i] for i in range(len(res))]+ [npv_mean]))

# print("fit flat PU template with nPV reweighting")
# yieldfitter.update_sigModel(sigTemplates_DY_flatPU+ "template_ZYield.root", hPV)
# res = yieldfitter.fit_simultanious(hZReco, hZReco_ss, "H_1")
# result.append(([res[i] for i in range(len(res))]+ [npv_mean]))
#
# print("fit high PU template with non Isolated Trigger w/o nPV reweighting")
# yieldfitter.update_sigModel(sigTemplates_DY_L1SingleMu25+ "template_ZYield.root")
# res = yieldfitter.fit_simultanious(hZReco, hZReco_ss, "H_2")
# result.append(([res[i] for i in range(len(res))]+ [npv_mean]))
#
# print("fit high PU template w/o nPV reweighting")
# yieldfitter.update_sigModel(sigTemplates_DY+ "template_ZYield.root")
# res = yieldfitter.fit_simultanious(hZReco, hZReco_ss, "H_3")
# result.append(([res[i] for i in range(len(res))]+ [npv_mean]))


results = {}
results["H"] = pd.DataFrame(result,
    columns=["fr","fr_err","tf","tf_err","chi2_OS","Chi2_SS","avgPU"])

hPV.Delete()
hZReco.Delete()
hZReco_ss.Delete()


eras = ["B", "C", "D", "E", "F"]
iStep = 120000
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

        # # new nPV distribution
        # hPV = TH1D("hPV_data", "th1d number of primary vertices shape", 75, -0.5, 74.5)
        # for nPV in selected['nPV']:
        #     hPV.Fill(nPV)
        # yieldfitter.update_sigModel(sigTemplates_DY+ "template_ZYield.root", hPV)
        # hPV.Delete()

        res = yieldfitter.fit_simultanious(hZReco_os_highPu, hZReco_ss_highPu, "{0}_{1}".format(era,iMeasure))

        result.append(([res[i] for i in range(len(res))] + [npv_mean]))

        iMeasure += 1
        hZReco_os_highPu.Delete()
        hZReco_ss_highPu.Delete()

    results[era] = pd.DataFrame(result, columns=["fr","fr_err","tf","tf_err","chi2_OS","Chi2_SS","avgPV"])


import pickle
pickle.dump(results, file(output+"/results.p","w"))
