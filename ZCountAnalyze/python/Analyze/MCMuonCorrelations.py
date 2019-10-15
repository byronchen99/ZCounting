from __future__ import division, print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from root_numpy import root2array, list_trees
from Utils import tree_to_df, plot_scatter
from scipy.stats import pearsonr
import pdb

import argparse

# plt.rcParams.update({'font.size': 18})

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
selection = 'z_genMass > 66 ' \
            '& z_genMass < 116 ' \
            '& muon_genPt > 27 ' \
            '& antiMuon_genPt > 27 ' \
            '& abs(muon_genEta) < 2.4 ' \
            '& abs(antiMuon_genEta) < 2.4 ' \
            '& muon_recoMatches == 1' \
            '& antiMuon_recoMatches == 1'

# specify which branches to load
branches = ['nPU',
            'muon_genEta', 'antiMuon_genEta',
            'muon_genPhi', 'antiMuon_genPhi',
            'muon_dxy', 'antiMuon_dxy',
            'muon_dz', 'antiMuon_dz',
            'muon_tkIso', 'antiMuon_tkIso',
            'muon_pfIso', 'antiMuon_pfIso'
            ]

print(">>> Load Events in gen acceptance")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
df = pd.concat(df)

print(">>> add new columns")
df['delRLL'] = np.sqrt(
    (df['muon_genEta'] - df['antiMuon_genEta']) ** 2 + (df['muon_genPhi'] - df['antiMuon_genPhi']) ** 2)


df = df.query('delRLL > 0.4')

plot_scatter(df['muon_tkIso'], df['antiMuon_tkIso'], '$\mu^{-}\ Tracker\ Isolation$', '$\mu^{+}\ Tracker\ Isolation$',
             range=(0, 1.), title='CMS Simulation',
             saveas='MuMu_inclusive_tkIso.png')
plot_scatter(df['muon_pfIso'], df['antiMuon_pfIso'], '$\mu^{-}\ PF\ Isolation$', '$\mu^{+}\ PF\ Isolation$',
             range=(0, 1.), title='CMS Simulation',
             saveas='MuMu_inclusive_pfIso.png')

exit()

plot_scatter(df['muon_dxy'], df['antiMuon_dxy'], '$\mu^{-} d_{xy}$', '$\mu^{+} d_{xy}$',
             range=(0, 0.5), title='CMS Simulation',
             saveas='MuMu_inclusive_dxy.png')
plot_scatter(df['muon_dz'], df['antiMuon_dz'], '$\mu^{-} d_{z}$', '$\mu^{+} d_{z}$',
             range=(0, 1.), title='CMS Simulation',
             saveas='MuMu_inclusive_dz.png')

dfLPU = df.query('nPU < 30')

plot_scatter(dfLPU['muon_dxy'], dfLPU['antiMuon_dxy'], '$\mu^{-} d_{xy}$', '$\mu^{+} d_{xy}$',
             range=(0, 0.5), cutsAdditional='nPU < 30',
             saveas='MuMu_LowPU_dxy.png')
plot_scatter(dfLPU['muon_dz'], dfLPU['antiMuon_dz'], '$\mu^{-} d_{z}$', '$\mu^{+} d_{z}$',
             range=(0, 1.), cutsAdditional='nPU < 30',
             saveas='MuMu_LowPU_dz.png')

dfHPU = df.query('nPU > 30')

plot_scatter(dfHPU['muon_dxy'], dfHPU['antiMuon_dxy'], '$\mu^{-} d_{xy}$', '$\mu^{+} d_{xy}$',
             range=(0, 0.5), cutsAdditional='nPU > 30',
             saveas='MuMu_HighPU_dxy.png')
plot_scatter(dfHPU['muon_dz'], dfHPU['antiMuon_dz'], '$\mu^{-} d_{z}$', '$\mu^{+} d_{z}$',
             range=(0, 1.), cutsAdditional='nPU > 30',
             saveas='MuMu_HighPU_dz.png')

dfBB = df.query('muon_genEta < 0.9 & antiMuon_genEta < 0.9')
dfEE = df.query('muon_genEta > 0.9 & antiMuon_genEta > 0.9')

plot_scatter(dfBB['muon_dxy'], dfBB['antiMuon_dxy'], '$\mu^{-} d_{xy}$', '$\mu^{+} d_{xy}$',
             range=(0, 0.5),
             cutsPtEta='p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 0.9',
             saveas='MuMu_BB_dxy.png')
plot_scatter(dfBB['muon_dz'], dfBB['antiMuon_dz'], '$\mu^{-} d_{z}$', '$\mu^{+} d_{z}$',
             range=(0, 1.),
             cutsPtEta='p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 0.9',
             saveas='MuMu_BB_dz.png')

plot_scatter(dfEE['muon_dxy'], dfEE['antiMuon_dxy'], '$\mu^{-} d_{xy}$', '$\mu^{+} d_{xy}$',
             range=(0, 0.5),
             cutsPtEta='p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad 0.9 < |\eta(\mu)| < 2.4',
             saveas='MuMu_EE_dxy.png')
plot_scatter(dfEE['muon_dz'], dfEE['antiMuon_dz'], '$\mu^{-} d_{z}$', '$\mu^{+} d_{z}$',
             range=(0, 1.),
             cutsPtEta='p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad 0.9 < |\eta(\mu)| < 2.4',
             saveas='MuMu_EE_dz.png')
