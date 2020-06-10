from __future__ import division, print_function

import numpy as np
import pandas as pd
import os
from root_numpy import root2array, list_trees
import matplotlib.pyplot as plt
import pdb

#local imports
from Utils import plot_scatter
os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df

import argparse

# Use TnPPairTrees

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
selection = 'is2HLT == 1 ' \
            #'| isSel == 1'

# specify which branches to load
branches = ['nPV','eventNumber',
            'eta1', 'eta2',
            'phi1', 'phi2',
            'dxy1', 'dxy2',
            'dz1', 'dz2',
            'pfIso1', 'pfIso2',
            'tkIso1', 'tkIso2',
            ]

print(">>> Load Events")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
df = pd.concat(df)

print(">>> add new columns")
df['delRLL'] = np.sqrt(
    (df['eta1'] - df['eta2']) ** 2 + (df['phi1'] - df['phi2']) ** 2)

print(">>> make delR cut")
df = df.query('delRLL > 0.4')

print("Number of events = {0}".format(len(df)))


def correlations(cut1, cut2):
    print("=====================================================")
    print("efficiency under cuts: {0} and {1}".format(cut1, cut2))
    effTag = len(df.query('{0}'.format(cut1)))/len(df)
    effProbe = len(df.query('{0}'.format(cut2)))/len(df)
    effProbe_Tag = len(df.query('{0} & {1}'.format(cut1, cut2)))/len(df.query('{0}'.format(cut1)))

    print("eff(Tag) = {0}".format(effTag))
    print("eff(Probe) = {0}".format(effProbe))
    print("eff(Probe|Tag) = {0}".format(effProbe_Tag))
    print("eff(Tag) * eff(Probe) = {0}".format(effTag * effProbe))
    print("eff(Tag) * eff(Probe|Tag) = {0}".format(effTag * effProbe_Tag))
    print("diff = {0} %".format(100* (effTag * effProbe_Tag - effTag * effProbe)))


correlations("dz1 < 0.5", "dz2 < 0.5")
correlations("dxy1 < 0.2", "dxy2 < 0.2")
correlations("tkIso1 < 0.1", "tkIso2 < 0.1")
correlations("pfIso1 < 0.2", "pfIso2 < 0.2")


plt.hist(df.drop_duplicates(subset='eventNumber')['nPV'], bins=np.linspace(0.5, 12.5, 13), histtype='step', color='k')
plt.xlabel('nPV')
plt.ylabel(r'#Events')
plt.savefig(output+"/nPV.png")
plt.clf()

print(">>> plot")
plot_scatter(df['tkIso1'], df['tkIso2'], '$\mu^{\mathrm{tag}}\ \mathrm{Tracker\ Isolation}$', '$\mu^{\mathrm{probe}}\ \mathrm{Tracker\ Isolation}$',
             range=(0, 1.), title='2017H (lowPU)',
             saveas=output+'/MuMu_inclusive_tkIso.png')
plot_scatter(df['pfIso1'], df['pfIso2'], '$\mu^{\mathrm{tag}}\ \mathrm{PF\ Isolation}$', '$\mu^{\mathrm{probe}}\ \mathrm{PF\ Isolation}$',
             range=(0, 1.), title=r'2017H (lowPU)',
             saveas=output+'/MuMu_inclusive_pfIso.png')

plot_scatter(df['dxy1'], df['dxy2'], '$\mu^{\mathrm{tag}} d_{xy}$', '$\mu^{\mathrm{probe}} d_{xy}$',
             range=(0, 0.5), title='2017H (lowPU)',
             saveas=output+'/MuMu_inclusive_dxy.png')
plot_scatter(df['dz1'], df['dz2'], '$\mu^{\mathrm{tag}} d_{z}$', '$\mu^{\mathrm{probe}} d_{z}$',
             range=(0, 1.), title='2017H (lowPU)',
             saveas=output+'/MuMu_inclusive_dz.png')
