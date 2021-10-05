from __future__ import division, print_function
import pandas as pd
import os
import numpy as np
from root_numpy import root2array, list_trees, array2tree
import pdb
import math

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
    '-y', '--year', default="2017", type=str,
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
            '& antiMuon_genRecoMatches == 1 ' \
            # '& (decayMode == 13 | decayMode == 151313)'


# specify which branches to load
branches = ['nPV', 'nPU',
            'eventweight',
            'Muon_eta[muon_genRecoObj]',
            'Muon_phi[muon_genRecoObj]',
            'Muon_pt[muon_genRecoObj]',
            'Muon_ID[muon_genRecoObj]',
            'Muon_isStandalone[muon_genRecoObj]',
            'Muon_nPixelHits[muon_genRecoObj]',
            'Muon_nTrackerLayers[muon_genRecoObj]',
            # 'Muon_isTracker[muon_genRecoObj]',
            'Muon_isGlobal[muon_genRecoObj]',
            'Muon_triggerBits[muon_genRecoObj]',
            # 'Muon_tkRelIso[muon_genRecoObj]',
            # 'Muon_pfRelIso04_all[muon_genRecoObj]',
            'Muon_eta[antiMuon_genRecoObj]',
            'Muon_phi[antiMuon_genRecoObj]',
            'Muon_pt[antiMuon_genRecoObj]',
            'Muon_ID[antiMuon_genRecoObj]',
            'Muon_isStandalone[antiMuon_genRecoObj]',
            'Muon_nPixelHits[antiMuon_genRecoObj]',
            'Muon_nTrackerLayers[antiMuon_genRecoObj]',
            # 'Muon_isTracker[antiMuon_genRecoObj]',
            'Muon_isGlobal[antiMuon_genRecoObj]',
            'Muon_triggerBits[antiMuon_genRecoObj]',
            # 'Muon_tkRelIso[antiMuon_genRecoObj]',
            # 'Muon_pfRelIso04_all[antiMuon_genRecoObj]',
            'z_recoMass',
            'decayMode'
            ]

print(">>> Load Events in acceptance")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 1) for i in inputs]
df = pd.concat(df)

df = df.rename(columns={
    'Muon_eta[muon_genRecoObj]_0':     'muon_recoEta',
    'Muon_eta[antiMuon_genRecoObj]_0': 'antiMuon_recoEta',
    'Muon_phi[muon_genRecoObj]_0':     'muon_recoPhi',
    'Muon_phi[antiMuon_genRecoObj]_0': 'antiMuon_recoPhi',
    'Muon_pt[muon_genRecoObj]_0':      'muon_recoPt',
    'Muon_pt[antiMuon_genRecoObj]_0':  'antiMuon_recoPt',
    'Muon_isStandalone[muon_genRecoObj]_0':     'muon_isStandalone',
    'Muon_isStandalone[antiMuon_genRecoObj]_0': 'antiMuon_isStandalone',
    'Muon_isTracker[muon_genRecoObj]_0':        'muon_isTracker',
    'Muon_isTracker[antiMuon_genRecoObj]_0':    'antiMuon_isTracker',
    'Muon_isGlobal[muon_genRecoObj]_0':         'muon_isGlobal',
    'Muon_isGlobal[antiMuon_genRecoObj]_0':     'antiMuon_isGlobal',
    'Muon_ID[muon_genRecoObj]_0':               'muon_ID',
    'Muon_ID[antiMuon_genRecoObj]_0':           'antiMuon_ID',
    'Muon_triggerBits[muon_genRecoObj]_0':      'muon_triggerBits',
    'Muon_triggerBits[antiMuon_genRecoObj]_0':  'antiMuon_triggerBits',
    'Muon_nPixelHits[muon_genRecoObj]_0':           'muon_nPixelHits',
    'Muon_nPixelHits[antiMuon_genRecoObj]_0':       'antiMuon_nPixelHits',
    'Muon_nTrackerLayers[muon_genRecoObj]_0':       'muon_nTrackerLayers',
    'Muon_nTrackerLayers[antiMuon_genRecoObj]_0':   'antiMuon_nTrackerLayers',
    # 'Muon_tkRelIso[muon_genRecoObj]_0': 'muon_tkIso',
    # 'Muon_tkRelIso[antiMuon_genRecoObj]_0': 'antiMuon_tkIso',
    # 'Muon_pfRelIso04_all[muon_genRecoObj]_0': 'muon_pfIso',
    # 'Muon_pfRelIso04_all[antiMuon_genRecoObj]_0': 'antiMuon_pfIso',
    'z_recoMass':'mass'
    })

df = df[abs(df['muon_recoEta']) < 2.4]
df = df[abs(df['antiMuon_recoEta']) < 2.4]
df = df[df['muon_recoPt'] > 27]
df = df[df['antiMuon_recoPt'] > 27]

df = df.astype({'mass': np.float64})

print(">>> add new columns")
df['delPhi'] = abs(
   abs(df['muon_recoPhi'] - df['antiMuon_recoPhi']).apply(lambda x: x - 2 * math.pi if x >= math.pi else x))
df['delR'] = np.sqrt(
   (df['muon_recoEta'] - df['antiMuon_recoEta']) ** 2 + df['delPhi'] ** 2)

df = df.query('delR > 0.8')

tkIsoCut = None #0.05
pfIsoCut = None #0.12
muonID = 4
#   1: Trk
#   2: Sta
#   3: Glo
#   4: tight ID w/o dxy and dz cuts
#   5: tight ID

print(">>> apply triggerfilter")
if args.year == "2016":
    triggerfilter = lambda x: 1 if (x&4 != 0) or (x&8 != 0) else 0    # Iso(Tk)Mu24

elif args.year == "2017":
    triggerfilter = lambda x: 1 if x&16 != 0 else 0  # IsoMu27

elif args.year == "2017H":
    triggerfilter = lambda x: 1 if x&256 != 0 else 0  # Mu17

elif args.year == "2018":
    triggerfilter = lambda x: 1 if x&4 != 0 else 0   # IsoMu24

else:
    print("ERROR: unknown year for triggerfilter")
    triggerfilter = lambda x: 0

df['muon_hlt'] = df['muon_triggerBits'].apply(triggerfilter)
df['antiMuon_hlt'] = df['antiMuon_triggerBits'].apply(triggerfilter)

print(">>> select events with a tag")
sel1 = 'muon_hlt == 1 & muon_ID >= {0}'.format(muonID)
sel2 = 'antiMuon_hlt == 1 & antiMuon_ID >= {0}'.format(muonID)
if tkIsoCut:
    sel1 += '& muon_tkIso < {0}'.format(tkIsoCut)
    sel2 += '& antiMuon_tkIso < {0}'.format(tkIsoCut)
if pfIsoCut:
    sel1 += '& muon_pfIso < {0}'.format(pfIsoCut)
    sel2 += '& antiMuon_pfIso < {0}'.format(pfIsoCut)

df = df.query('( {0} ) | ( {1} )'.format(sel1, sel2))
# <<<


# Event weight manipulation

print(">>> Muon selection")
df['passSel'] = ((df['muon_ID'] >= muonID) & (df['antiMuon_ID'] >= muonID))

if tkIsoCut:
    print(">>> degrade selected muons for isolation")
    df['passSel'] = (df['passSel']
        & (df['muon_tkIso'] < tkIsoCut) & (df['antiMuon_tkIso'] < tkIsoCut))

if pfIsoCut:
    print(">>> degrade selected muons for isolation")
    df['passSel'] = (df['passSel']
        & (df['muon_pfIso'] < pfIsoCut) & (df['antiMuon_pfIso'] < pfIsoCut))

print(">>> HLT selection")
df['passHLT'] = (df['passSel']
    & (df['muon_hlt'] == 1) & (df['antiMuon_hlt'] == 1))

df['isSta'] = ((df['muon_isStandalone']) & (df['antiMuon_isStandalone']))
# df['isTrk']   = ((df['muon_isTracker'])    & (df['antiMuon_isTracker']))
df['isTrk']   = ((df['muon_nPixelHits'] >= 1)    & (df['antiMuon_nPixelHits'] >= 1)
    & (df['muon_nTrackerLayers'] >= 6)    & (df['antiMuon_nTrackerLayers'] >= 6))

# df['isTrk'] = ((df['muon_nPixelHits'] >= 0) & (df['antiMuon_nPixelHits'] >= 0)
#     & (df['muon_nTrackerLayers'] >= 0) & (df['antiMuon_nTrackerLayers'] >= 0))

df['passGlo'] = ((df['isTrk']) & (df['muon_isGlobal']) & (df['antiMuon_isGlobal']))


print(">>> Template for 2HLT and 1HLT")

df2HLT = df.query('passHLT == 1')
df2HLT = df2HLT.rename(columns={'muon_recoPt': 'ptTag',
                                'antiMuon_recoPt': 'ptProbe',
                                'muon_recoEta': 'etaTag',
                                'antiMuon_recoEta': 'etaProbe'
                                })
df1HLT_1 = df.query('passHLT == 0 & passSel == 1 & muon_hlt == 0')
df1HLT_1 = df1HLT_1.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
df1HLT_2 = df.query('passHLT == 0 & passSel == 1 & antiMuon_hlt == 0')
df1HLT_2 = df1HLT_2.rename(columns={'muon_recoPt': 'ptTag',
                                        'antiMuon_recoPt': 'ptProbe',
                                        'muon_recoEta': 'etaTag',
                                        'antiMuon_recoEta': 'etaProbe'
                                        })

dfOut = pd.concat([df2HLT, df1HLT_1, df1HLT_2], sort=True)
dfOut = dfOut.rename(columns={'z_recoMass': 'mass', 'passHLT': 'pass'})
dfOut = dfOut.astype({'pass': np.bool})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass','eventweight']]

tfile = ROOT.TFile.Open(output+"/template_ZYield.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
hPV.Write()

ttree.Delete()

tfile.Close()

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
dfOut = dfOut.astype({'pass': np.bool})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass','eventweight']]

tfile = ROOT.TFile.Open(output+"/template_Sel.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
hPV.Write()
ttree.Delete()

tfile.Close()

### --- global to track -> fail outer track (standalone)
print(">>> Template for standalone")
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
dfOut = dfOut.astype({'pass': np.bool})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass','eventweight']]

tfile = ROOT.TFile.Open(output+"/template_Sta.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
ttree.Delete()
hPV.Write()

tfile.Close()

### --- global to standalone -> fail inner track
print(">>> Template for inner track")
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
dfOut = dfOut.astype({'pass': np.bool})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass','eventweight']]

tfile = ROOT.TFile.Open(output+"/template_Trk.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
ttree.Delete()
hPV.Write()

tfile.Close()

hPV.Delete()
