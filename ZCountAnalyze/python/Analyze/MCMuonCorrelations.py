from __future__ import division, print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from root_numpy import root2array, list_trees
from Utils import plot_scatter
from scipy.stats import pearsonr
import pdb
import math

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df

import argparse

# plt.rcParams.update({'font.size': 18})

parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-i', '--input', type=str,
    help='specify input hdf5 file'
)
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output dir'
)
args = parser.parse_args()

input = args.input
output = args.output

print(">>> load dataframe")
store = pd.HDFStore(input)
df = store['dfGen']  # load it

if not os.path.isdir(output):
    os.mkdir(output)

# print(">>> add new columns")
# df['delPhiLL'] = abs(df['muon_genPhi'] - df['antiMuon_genPhi']).apply(lambda x: x - 2 * math.pi if x > math.pi else x)

    # (df['muon_genEta'] - df['antiMuon_genEta']) ** 2 + (df['muon_genPhi'] - df['antiMuon_genPhi']) ** 2)
# df['delEtaLL'] = abs(df['muon_genEta'] - df['antiMuon_genEta'])

# df['delRLL'] = np.sqrt( df['delPhiLL']**2 + df['delEtaLL']**2)
#
# df = df.query('delRLL > 0.4')

# dfFromPV = df.query('muon_isFromPV == 1 & antiMuon_isFromPV == 1')
# dfNotFromPV = df.query('muon_isFromPV == 0 & antiMuon_isFromPV == 0')
#
# plot_scatter(dfFromPV['muon_dxy'], dfFromPV['antiMuon_dxy'], '$\mu^{-} d_{xy}$', '$\mu^{+} d_{xy}$',
#              range=(0, 0.5), title='Events with Z from PV',
#              saveas=output+'/MuMu_FromPV_dxy.png')
# plot_scatter(dfFromPV['muon_dz'], dfFromPV['antiMuon_dz'], '$\mu^{-} d_{z}$', '$\mu^{+} d_{z}$',
#              range=(0, .5), title='Events with Z from PV',
#              saveas=output+'/MuMu_FromPV_dz.png')
#
# plot_scatter(dfNotFromPV['muon_dxy'], dfNotFromPV['antiMuon_dxy'], '$\mu^{-} d_{xy}$', '$\mu^{+} d_{xy}$',
#              range=(0, 0.5), title='Events with Z not from PV',
#              saveas=output+'/MuMu_NotFromPV_dxy.png')
# plot_scatter(dfNotFromPV['muon_dz'], dfNotFromPV['antiMuon_dz'], '$\mu^{-} d_{z}$', '$\mu^{+} d_{z}$',
#              range=(0, .5), title='Events with Z not from PV',
#              saveas=output+'/MuMu_NotFromPV_dz.png')
#
# exit()
# pfIso = 0.12
# tkIso = 0.05
#
# df = df.loc[df['muon_ID'] >= 4]
# df = df.loc[df['antiMuon_ID'] >= 4]

# df = df.loc[df['muon_pfIso'] < pfIso]
# df = df.loc[df['antiMuon_pfIso'] < pfIso]
#
# df = df.loc[df['muon_tkIso'] < tkIso]
# df = df.loc[df['antiMuon_tkIso'] < tkIso]


# plot_scatter(df['muon_tkIso'], df['antiMuon_tkIso'], r'$\mu^{-}$ Tracker isolation', '$\mu^{+}$ Tracker isolation',
#              range=(0, 3.), #title='CMS Simulation',
#              saveas=output+'/MuMu_inclusive_tkIso.png')
# plot_scatter(df['muon_pfIso'], df['antiMuon_pfIso'], r'$\mu^{-}$ Particle flow isolation', '$\mu^{+}$ Particle flow isolation',
#              range=(0, 3.), #title='CMS Simulation',
#              saveas=output+'/MuMu_inclusive_pfIso.png')
#
# plot_scatter(df['muon_dxy'], df['antiMuon_dxy'], '$\mu^{-} d_\mathrm{XY}$ [cm]', '$\mu^{+} d_\mathrm{XY}$ [cm]',
#              range=(0, 3.0), #title='CMS Simulation',
#              saveas=output+'/MuMu_inclusive_dxy.png')
# plot_scatter(df['muon_dz'], df['antiMuon_dz'], '$\mu^{-} d_\mathrm{Z}$ [cm]', '$\mu^{+} d_\mathrm{Z}$ [cm]',
#              range=(0, 3.0), #title='CMS Simulation',
#              saveas=output+'/MuMu_inclusive_dz.png')

plot_scatter(df['muon_pt'], df['antiMuon_pt'], '$\mu^{-} p_\mathrm{T}$ [GeV]', '$\mu^{+} p_\mathrm{T}$ [GeV]',
             range=(15, 80), eventWeights=df['eventweight'],
             #title='CMS Simulation',
             saveas=output+'/MuMu_inclusive_pt.pdf')
plot_scatter(df['muon_phi'], df['antiMuon_phi'], '$\mu^{-} \phi$', '$\mu^{+} \phi$',
             range=(-math.pi, math.pi), eventWeights=df['eventweight'], #title='CMS Simulation',
             saveas=output+'/MuMu_inclusive_phi.pdf')
plot_scatter(df['muon_eta'], df['antiMuon_eta'], '$\mu^{-} \eta$', '$\mu^{+} \eta$',
             range=(-2.4, 2.4), eventWeights=df['eventweight'], #title='CMS Simulation',
             saveas=output+'/MuMu_inclusive_eta.pdf')
# plot_scatter(df['delPhiLL'], df['delEtaLL'], '$\Delta \phi (\mu^{-},\mu^{+})$', '$\Delta \eta (\mu^{-},\mu^{+})$',
#              range=(0, math.pi), #title='CMS Simulation',
#              rangey=(0,5.),
#              saveas=output+'/MuMu_inclusive_delPhiEta.pdf')
# plot_scatter(df['delRLL'], df['z_genMass'], '$m(\mu^{-},\mu^{+})$', '$\Delta R(\mu^{-},\mu^{+})$',
#              range=(86, 96), #title='CMS Simulation',
#              rangey=(0,5.),
#              saveas=output+'/MuMu_inclusive_delRMass.pdf')

exit()
# >>> now reqire HLT_IsoMu27 at one of the muons (mu^minus)

df = df.loc[df['muon_hlt_4'] == 1]

plot_scatter(df['muon_tkIso'], df['antiMuon_tkIso'], '$\mu^{-}\ Tracker\ Isolation$', '$\mu^{+}\ Tracker\ Isolation$',
             range=(0, 3.), #title='CMS Simulation',
             saveas=output+'/MuMu_inclusive_tkIso_MuHLTIsoMu27.png')
plot_scatter(df['muon_pfIso'], df['antiMuon_pfIso'], '$\mu^{-}\ PF\ Isolation$', '$\mu^{+}\ PF\ Isolation$',
             range=(0, 3.), #title='CMS Simulation',
             saveas=output+'/MuMu_inclusive_pfIso_MuHLTIsoMu27.png')


plot_scatter(df['muon_dxy'], df['antiMuon_dxy'], '$\mu^{-} d_{xy}$', '$\mu^{+} d_{xy}$',
             range=(0, 3.0), #title='CMS Simulation',
             saveas=output+'/MuMu_inclusive_dxy_MuHLTIsoMu27.png')
plot_scatter(df['muon_dz'], df['antiMuon_dz'], '$\mu^{-} d_{z}$ [cm]', '$\mu^{+} d_{z}$ [cm]',
             range=(0, 3.0), #title='CMS Simulation',
             saveas=output+'/MuMu_inclusive_dz_MuHLTIsoMu27.png')


exit()

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
