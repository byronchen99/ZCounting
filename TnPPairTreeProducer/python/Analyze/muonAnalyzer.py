from __future__ import division, print_function

# production of histograms from TnP trees
#   the histograms have the same format as the onse produced by DQM
#   and can therefore used with the same tools, e.g. ZCounting.py

import numpy as np
import pandas as pd
import os
import root_numpy as rn
import pdb
import ROOT
import argparse
import glob
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
#local imports
os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df


parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output dir'
)

args = parser.parse_args()

output = args.output


if not os.path.isdir(output):
    os.mkdir(output)

####################
# loading of input #
####################

input = glob.glob("/pnfs/desy.de/cms/tier2/store/user/dwalter/SingleMuon/MuonTrees_V06/200610_143117/0000/output_0_*")

treeName = rn.list_trees(input[0])[0]
print(">>> Load Events from {0} files".format(len(input)))
_df = []
for i, ifile in enumerate(input):
    print("> file {0}/{1}".format(i, len(input)))
    sys.stdout.flush()
    tfile = ROOT.TFile(ifile)
    ttree = tfile.Get(treeName)
    if ttree.GetEntries() == 0:
        print("tree has no entries, continue")
        continue
    _df.append(tree_to_df(rn.root2array(ifile, treeName)))

print(">>> Concatenate")
df = pd.concat(_df)

###############################################
# computation of efficiencies and acceptances #
###############################################

def effacc(selection):
    """
    selection:
        lumbda function that returns a bool like lambda x: df['pt'] > x
    return:
        function for efficiency depending on selection as defined in lambda function
        eff(x)
    """
    _eff = lambda x: len(df.loc[(df['is_HLT_IsoMu27'] == 1) & (selection(x))]) / len(df.loc[selection(x)])
    _acc = lambda x: len(df.loc[(df['is_HLT_IsoMu27'] == 1) & (selection(x))]) / len(df.loc[df['is_HLT_IsoMu27'] == 1])

    return np.vectorize(_eff), np.vectorize(_acc)

def effacc2d(selection_x, selection_y):
    """
    selection_{x,y}:
        lumbda function that returns a bool like lambda x: df['pt'] > x
    return:
        function for efficiency depending on selections as defined in lambda function
        eff(x, y)
    """
    _eff = lambda x,y: len(df.loc[(df['is_HLT_IsoMu27'] == 1) & (selection_x(x)) & (selection_y(y))]) / len(df.loc[selection_x(x) & (selection_y(y))])
    _acc = lambda x,y: len(df.loc[(df['is_HLT_IsoMu27'] == 1) & (selection_x(x)) & (selection_y(y))]) / len(df.loc[df['is_HLT_IsoMu27'] == 1])

    return np.vectorize(_eff), np.vectorize(_acc)

eff_pt, acc_pt = effacc(lambda x: df['pt'] > x)
#eff_eta, acc_eta = effacc(lambda x: df['eta'] < x)

#eff_pfIso, acc_pfIso = effacc(lambda x: df['pfIso'] < x)
#eff_tkIso, acc_tkIso = effacc(lambda x: df['tkIso'] < x)

eff_tkIso_pfIso, acc_tkIso_pfIso = effacc2d(lambda x: df['tkIso'] < x, lambda x: df['pfIso'] < x)
#eff_dxy_dz, acc_dxy_dz = effacc2d(lambda x: df['dxy'] < x, lambda x: df['dz'] < x)

############
# plotting #
############

print(">>> plot")

def plot(range, eff, acc, xlabel='', plot_max_eff=False):
    xx = np.linspace(range[0], range[1], 200)
    plt.cla()
    yy_eff = eff(xx)
    yy_acc = acc(xx)

    plt.plot(xx, yy_eff, label='efficiency')
    plt.plot(xx, yy_acc, label='acceptance')

    if plot_max_eff:
        eff_max = xx[np.argmax(yy_eff)]
        ymin = min(min(yy_eff), min(yy_acc))
        ymax = max(max(yy_eff), max(yy_acc))
        plt.plot([eff_max,eff_max], [ymin,ymax], '-k', label='max(efficiency)')

    plt.xlabel(xlabel)
    plt.legend()
    plt.savefig(output+"/triggerdiff_{0}_{1}to{2}.png".format(xlabel, str(range[0]).replace(".",""), str(range[1]).replace(".","")))

def plot2d(xrange, yrange, eff, acc, xlabel='', ylabel='', suffix=''):

    def _plot2d(func, zlabel):
        xx, yy = np.meshgrid(np.linspace(xrange[0], xrange[1], 20), np.linspace(yrange[0], yrange[1], 20))#, sparse=True)

        plt.cla()
        fig, ax = plt.subplots()
        cs = ax.contourf(xx, yy, func(xx, yy), cmap=cm.PuBu_r)
        cbar = fig.colorbar(cs)
        cbar.ax.set_ylabel(zlabel)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(suffix)
        plt.savefig(output+"/triggerdiff2d_{0}{1}to{2}_{3}{4}to{5}_{6}{7}.png".format(
            xlabel, str(xrange[0]).replace(".",""), str(xrange[1]).replace(".",""),
            ylabel, str(yrange[0]).replace(".",""), str(yrange[1]).replace(".",""),
            zlabel, suffix))

    _plot2d(eff, 'efficiency')
    _plot2d(acc, 'acceptance')


plot((27, 55), eff_pt, acc_pt, xlabel="pt")
plot((27, 30), eff_pt, acc_pt, xlabel="pt")

# plot((0.5, 2.4), eff_pt, acc_pt, xlabel="eta")
# plot((0.025, 0.2), eff_tkIso, acc_tkIso, xlabel="tkIso")
# plot((0.06, 0.4), eff_pfIso, acc_pfIso, xlabel="pfIso")
#
# plot((0.01, 1), eff_tkIso, acc_tkIso, xlabel="tkIso")
# plot((0.01, 1), eff_pfIso, acc_pfIso, xlabel="pfIso")

plot2d((0.025, 0.2), (0.06, 0.4), eff_tkIso_pfIso, acc_tkIso_pfIso, 'tkIso', 'pfIso')

df = df.query('pt > 28')
plot2d((0.025, 0.2), (0.06, 0.4), eff_tkIso_pfIso, acc_tkIso_pfIso, 'tkIso', 'pfIso','_pt28')

df = df.query('pt > 30')
plot2d((0.025, 0.2), (0.06, 0.4), eff_tkIso_pfIso, acc_tkIso_pfIso, 'tkIso', 'pfIso','_pt30')

df = df.query('pt > 32')
plot2d((0.025, 0.2), (0.06, 0.4), eff_tkIso_pfIso, acc_tkIso_pfIso, 'tkIso', 'pfIso','_pt32')

#plot2d((0.1, 1), (0.25, 2), eff_dxy_dz, acc_dxy_dz, 'dxy', 'dz')
