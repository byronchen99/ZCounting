from __future__ import division, print_function
import pandas as pd
import os
from root_numpy import root2array, list_trees, array2tree
from Utils import tree_to_df

import argparse
import ROOT

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


# acceptance selection
selection = 'z_recoMass > 66 ' \
            '& z_recoMass < 116 ' \
            '& muon_recoPt > 27 ' \
            '& antiMuon_recoPt > 27 ' \
            '& abs(muon_recoEta) < 2.4 ' \
            '& abs(antiMuon_recoEta) < 2.4 ' \
            '& muon_recoMatches == 1 ' \
            '& antiMuon_recoMatches == 1'

# specify which branches to load
branches = ['nPV', 'nPU', 'z_recoMass',
            'muon_recoEta', 'antiMuon_recoEta',
            'muon_recoPhi', 'antiMuon_recoPhi',
            'muon_recoPt',  'antiMuon_recoPt',
            'muon_ID', 'antiMuon_ID',
            'muon_triggerBits', 'antiMuon_triggerBits',
            ]

print(">>> Load Events in acceptance")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
df = pd.concat(df)

#print(">>> add new columns")
#df['delRLL'] = np.sqrt(
#    (df['muon_genEta'] - df['antiMuon_genEta']) ** 2 + (df['muon_genPhi'] - df['antiMuon_genPhi']) ** 2)

#  df = df.query('delRLL > 0.4')

print(">>> convert bit code into bit map")
for iBit in range(0, 5):
    nBit = 2 ** iBit
    iBit += 1
    df['muon_hlt_{0}'.format(iBit)] = df['muon_triggerBits'].apply(
        lambda x: 1 if x % (nBit * 2) >= nBit else 0)
    df['antiMuon_hlt_{0}'.format(iBit)] = df['antiMuon_triggerBits'].apply(
        lambda x: 1 if x % (nBit * 2) >= nBit else 0)

print(">>> select events with a tag")

df = df.query('(muon_hlt_1 == 1 & muon_ID > 4) | (antiMuon_hlt_1 == 1 & antiMuon_ID > 4)')
df['passHLT'] = ((df['muon_ID'] >= 4) & (df['antiMuon_ID'] >= 4) & (df['muon_hlt_1'] == 1) & (df['antiMuon_hlt_1'] == 1))
df['passSel'] = ((df['muon_ID'] >= 4) & (df['antiMuon_ID'] >= 4))
df['passGlo'] = ((df['muon_ID'] >= 3) & (df['antiMuon_ID'] >= 3))
df['passTrk'] = ((df['muon_ID'] >= 1) & (df['antiMuon_ID'] >= 1))

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
dfFailHLT1 = df.query('passHLT == 0 & passSel == 1 & muon_hlt_1 == 0')
dfFailHLT1 = dfFailHLT1.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
dfFailHLT2 = df.query('passHLT == 0 & passSel == 1 & antiMuon_hlt_1 == 0')
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
ttree.Delete()

tfile.Close()


print(">>> Template for selection step")
dfFailSel1 = df.query('passSel == 0 & passGlo == 1 & muon_hlt_1 == 0')
dfFailSel1 = dfFailSel1.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
dfFailSel2 = df.query('passSel == 0 & passGlo == 1 & antiMuon_hlt_1 == 0')
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
ttree.Delete()

tfile.Close()

print(">>> Template for global step")
dfFailGlo1 = df.query('passGlo == 0 & passTrk == 1 & muon_hlt_1 == 0')
dfFailGlo1 = dfFailGlo1.rename(columns={'muon_recoPt': 'ptProbe',
                                        'antiMuon_recoPt': 'ptTag',
                                        'muon_recoEta': 'etaProbe',
                                        'antiMuon_recoEta': 'etaTag'
                                        })
dfFailGlo2 = df.query('passGlo == 0 & passTrk == 1 & antiMuon_hlt_1 == 0')
dfFailGlo2 = dfFailGlo2.rename(columns={'muon_recoPt': 'ptTag',
                                        'antiMuon_recoPt': 'ptProbe',
                                        'muon_recoEta': 'etaTag',
                                        'antiMuon_recoEta': 'etaProbe'
                                        })

dfOut = pd.concat([dfPassHLT1, dfPassHLT2, dfFailHLT1, dfFailHLT2, dfFailSel1, dfFailSel2, dfFailGlo1, dfFailGlo2], sort=True)
dfOut = dfOut.rename(columns={'z_recoMass': 'mass', 'passGlo': 'pass'})
dfOut = dfOut[['nPV', 'nPU', 'ptTag', 'ptProbe', 'etaTag', 'etaProbe', 'mass', 'pass']]

tfile = ROOT.TFile.Open(output+"/template_Glo.root", "RECREATE")

ttree = array2tree(dfOut.to_records(index=False))
ttree.Write()
ttree.Delete()

tfile.Close()
