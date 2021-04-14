from __future__ import division, print_function
import pandas as pd
import os
import numpy as np
from root_numpy import root2array, list_trees, array2tree
import pdb

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df

import argparse
import ROOT

parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-i', '--input', nargs='+',
    help='specify input root file'
)
parser.add_argument(
    '-y', '--year', default=2017, type=int,
    help='specify year'
)
parser.add_argument(
    '-o', '--output', type=str, default='./',
    help='specify output dir'
)
args = parser.parse_args()

inputs = args.input
output = args.output

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


dfPV = [tree_to_df(root2array(i, treeName[0], branches=['nPV']), 5) for i in inputs]
dfPV = pd.concat(dfPV)

hPV = ROOT.TH1D("hPV", "th1d number of primary vertices shape", 75, -0.5, 74.5)

for nPV in dfPV['nPV']:
    hPV.Fill(nPV)

hPV.Scale(1. / hPV.Integral())

# acceptance selection
selection = 'muon_genRecoMatches == 1 ' \
            '& antiMuon_genRecoMatches == 1 '


# specify which branches to load
branches = ['nPV', 'nPU',
            'Muon_eta[muon_genRecoObj]',
            'Muon_phi[muon_genRecoObj]',
            'Muon_pt[muon_genRecoObj]',
            'Muon_ID[muon_genRecoObj]',
            'Muon_triggerBits[muon_genRecoObj]',
            'Muon_tkRelIso[muon_genRecoObj]',
            'Muon_pfRelIso04_all[muon_genRecoObj]',
            'Muon_eta[antiMuon_genRecoObj]',
            'Muon_phi[antiMuon_genRecoObj]',
            'Muon_pt[antiMuon_genRecoObj]',
            'Muon_ID[antiMuon_genRecoObj]',
            'Muon_triggerBits[antiMuon_genRecoObj]',
            'Muon_tkRelIso[antiMuon_genRecoObj]',
            'Muon_pfRelIso04_all[antiMuon_genRecoObj]',
            'z_recoMass'
            ]

print(">>> Load Events in acceptance")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 1) for i in inputs]
df = pd.concat(df)

df = df.rename(columns={
    'Muon_eta[muon_genRecoObj]_0': 'muon_recoEta',
    'Muon_eta[antiMuon_genRecoObj]_0': 'antiMuon_recoEta',
    'Muon_phi[muon_genRecoObj]_0': 'muon_recoPhi',
    'Muon_phi[antiMuon_genRecoObj]_0': 'antiMuon_recoPhi',
    'Muon_pt[muon_genRecoObj]_0': 'muon_recoPt',
    'Muon_pt[antiMuon_genRecoObj]_0': 'antiMuon_recoPt',
    'Muon_ID[muon_genRecoObj]_0': 'muon_ID',
    'Muon_ID[antiMuon_genRecoObj]_0': 'antiMuon_ID',
    'Muon_triggerBits[muon_genRecoObj]_0': 'muon_triggerBits',
    'Muon_triggerBits[antiMuon_genRecoObj]_0': 'antiMuon_triggerBits',
    'Muon_tkRelIso[muon_genRecoObj]_0': 'muon_tkIso',
    'Muon_tkRelIso[antiMuon_genRecoObj]_0': 'antiMuon_tkIso',
    'Muon_pfRelIso04_all[muon_genRecoObj]_0': 'muon_pfIso',
    'Muon_pfRelIso04_all[antiMuon_genRecoObj]_0': 'antiMuon_pfIso',
    'z_recoMass':'mass'
    })

df = df[abs(df['muon_recoEta']) < 2.4]
df = df[abs(df['antiMuon_recoEta']) < 2.4]

df = df.astype({'mass': np.float64})

# print(">>> compute Z mass")
# df['mass'] = df.apply(lambda x: (
#     ROOT.Math.PtEtaPhiMVector(x['muon_recoPt'], x['muon_recoEta'], x['muon_recoPhi'], 0.105658369)
#     + ROOT.Math.PtEtaPhiMVector(x['antiMuon_recoPt'], x['antiMuon_recoEta'], x['antiMuon_recoPhi'], 0.105658369)).M(),
#     axis=1)

#print(">>> add new columns")
#df['delRLL'] = np.sqrt(
#    (df['muon_genEta'] - df['antiMuon_genEta']) ** 2 + (df['muon_genPhi'] - df['antiMuon_genPhi']) ** 2)

#  df = df.query('delRLL > 0.4')

tkIsoCut = 0.05#999999 #
pfIsoCut = 0.12
muonID = 5
#   1: Trk
#   2: Sta
#   3: Glo
#   4: tight ID w/o dxy and dz cuts
#   5: tight ID

print(">>> convert bit code into bit map")
if args.year == 2016:
    # For 2016:
    #   0: "HLT_IsoMu24_v*"
    #   1: "HLT_IsoTkMu24_v*"
    # -> Do an or of both triggers
    df['muon_hlt'] = df['muon_triggerBits'].apply(
        lambda x: 1 if x >= 1 else 0)
    df['antiMuon_hlt'] = df['antiMuon_triggerBits'].apply(
        lambda x: 1 if x >= 1 else 0)
else:
    # For 2017/2018:
    #   0: "HLT_L1SingleMu18_v*"
    #   1: "HLT_L1SingleMu25_v*"
    #   2: "HLT_IsoMu24_v*"
    #   3: "HLT_IsoMu27_v*"
    #   4: "HLT_IsoMu30_v*"
    iBit = 3

    nBit = 2 ** iBit
    df['muon_hlt'] = df['muon_triggerBits'].apply(
        lambda x: 1 if x % (nBit * 2) >= nBit else 0)
    df['antiMuon_hlt'] = df['antiMuon_triggerBits'].apply(
        lambda x: 1 if x % (nBit * 2) >= nBit else 0)

print(">>> select events with a tag")
df = df.query('(muon_hlt == 1 & muon_ID >= 4 & muon_tkIso < {0} & muon_pfIso < {1}) | \
    (antiMuon_hlt == 1 & antiMuon_ID >= 4 & antiMuon_tkIso < {0} & antiMuon_pfIso < {1})'.format(tkIsoCut, pfIsoCut))

df['passHLT'] = ((df['muon_ID'] >= muonID) & (df['antiMuon_ID'] >= muonID) \
    & (df['muon_hlt'] == 1) & (df['antiMuon_hlt'] == 1) \
    & (df['muon_tkIso'] < tkIsoCut) & (df['antiMuon_tkIso'] < tkIsoCut) \
    & (df['muon_pfIso'] < pfIsoCut) & (df['antiMuon_pfIso'] < pfIsoCut))

df['passSel'] = ((df['muon_ID'] >= muonID) & (df['antiMuon_ID'] >= muonID)\
    & (df['muon_tkIso'] < tkIsoCut) & (df['antiMuon_tkIso'] < tkIsoCut) \
    & (df['muon_pfIso'] < pfIsoCut) & (df['antiMuon_pfIso'] < pfIsoCut))

df['passGlo'] = ((df['muon_ID'] >= 3) & (df['antiMuon_ID'] >= 3))
df['isSta'] = ((df['muon_ID'] == 2) & (df['antiMuon_ID'] == 2))
df['isTrk'] = ((df['muon_ID'] == 1) & (df['antiMuon_ID'] == 1))

print(">>> Template for hlt step")
dfPassHLT1 = df.query('passHLT == 1')
dfPassHLT1 = dfPassHLT1.rename(columns={'muon_recoPt': 'ptTag',
                                        'antiMuon_recoPt': 'ptProbe',
                                        'muon_recoEta': 'etaTag',
                                        'antiMuon_recoEta': 'etaProbe'
                                        })
dfPassHLT2 = df.query('passHLT == 1')
dfPassHLT2 = dfPassHLT2.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
dfFailHLT1 = df.query('passHLT == 0 & passSel == 1 & muon_hlt == 0')
dfFailHLT1 = dfFailHLT1.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
dfFailHLT2 = df.query('passHLT == 0 & passSel == 1 & antiMuon_hlt == 0')
dfFailHLT2 = dfFailHLT2.rename(columns={'muon_recoPt': 'ptTag',
                                        'antiMuon_recoPt': 'ptProbe',
                                        'muon_recoEta': 'etaTag',
                                        'antiMuon_recoEta': 'etaProbe'
                                        })

dfOut = pd.concat([dfPassHLT1, dfPassHLT2, dfFailHLT1, dfFailHLT2], sort=True)
dfOut = dfOut.rename(columns={'z_recoMass': 'mass', 'passHLT': 'pass'})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass']]

tfile = ROOT.TFile.Open(output+"/template_HLT.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
hPV.Write()

ttree.Delete()

tfile.Close()


dfOut = pd.concat([dfPassHLT1, dfFailHLT1, dfFailHLT2], sort=True)
dfOut = dfOut.rename(columns={'z_recoMass': 'mass', 'passSel': 'pass'})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass']]

tfile = ROOT.TFile.Open(output+"/template_ZYield.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
hPV.Write()

ttree.Delete()

tfile.Close()


print(">>> Template for selection step")
dfFailSel1 = df.query('passSel == 0 & passGlo == 1 & muon_hlt == 0')
dfFailSel1 = dfFailSel1.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
dfFailSel2 = df.query('passSel == 0 & passGlo == 1 & antiMuon_hlt == 0')
dfFailSel2 = dfFailSel2.rename(columns={'muon_recoPt': 'ptTag',
                                        'antiMuon_recoPt': 'ptProbe',
                                        'muon_recoEta': 'etaTag',
                                        'antiMuon_recoEta': 'etaProbe'
                                        })

dfOut = pd.concat([dfPassHLT1, dfPassHLT2, dfFailHLT1, dfFailHLT2, dfFailSel1, dfFailSel2], sort=True)
dfOut = dfOut.rename(columns={'z_recoMass': 'mass', 'passSel': 'pass'})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass']]

tfile = ROOT.TFile.Open(output+"/template_Sel.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
hPV.Write()
ttree.Delete()

tfile.Close()

### --- global to track -> fail outer track (standalone)
print(">>> Template for global step")
dfFailSta1 = df.query('passGlo == 0 & isTrk == 1 & muon_hlt == 0')
dfFailSta1 = dfFailSta1.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
dfFailSta2 = df.query('passGlo == 0 & isTrk == 1 & antiMuon_hlt == 0')
dfFailSta2 = dfFailSta2.rename(columns={'muon_recoPt': 'ptTag',
                                        'antiMuon_recoPt': 'ptProbe',
                                        'muon_recoEta': 'etaTag',
                                        'antiMuon_recoEta': 'etaProbe'
                                        })

dfOut = pd.concat([dfPassHLT1, dfPassHLT2, dfFailHLT1, dfFailHLT2, dfFailSel1, dfFailSel2, dfFailSta1, dfFailSta2], sort=True)
dfOut = dfOut.rename(columns={'z_recoMass': 'mass', 'passGlo': 'pass'})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass']]

tfile = ROOT.TFile.Open(output+"/template_Sta.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
ttree.Delete()
hPV.Write()

tfile.Close()

### --- global to standalone -> fail inner track
print(">>> Template for global step")
dfFailTrk1 = df.query('passGlo == 0 & isSta == 1 & muon_hlt == 0')
dfFailTrk1 = dfFailTrk1.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
dfFailTrk2 = df.query('passGlo == 0 & isSta == 1 & antiMuon_hlt == 0')
dfFailTrk2 = dfFailTrk2.rename(columns={'muon_recoPt': 'ptTag',
                                        'antiMuon_recoPt': 'ptProbe',
                                        'muon_recoEta': 'etaTag',
                                        'antiMuon_recoEta': 'etaProbe'
                                        })

dfOut = pd.concat([dfPassHLT1, dfPassHLT2, dfFailHLT1, dfFailHLT2, dfFailSel1, dfFailSel2, dfFailTrk1, dfFailTrk2], sort=True)
dfOut = dfOut.rename(columns={'z_recoMass': 'mass', 'passGlo': 'pass'})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass']]

tfile = ROOT.TFile.Open(output+"/template_Trk.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
ttree.Delete()
hPV.Write()

tfile.Close()

hPV.Delete()
