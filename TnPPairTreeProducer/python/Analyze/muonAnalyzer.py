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


parser = argparse.ArgumentParser(prog='./difficiencies')
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
# computation of difficiencies and acceptances #
###############################################

def diffacc(selection):
    """
    selection:
        lumbda function that returns a bool like lambda x: df['pt'] > x
    return:
        function for difficiency depending on selection as defined in lambda function
        diff(x)
    """
    _diff = lambda x: len(df.loc[(df['is_HLT_IsoMu27'] == 1) & (selection(x))]) / len(df.loc[selection(x)])
    _acc = lambda x: len(df.loc[(df['is_HLT_IsoMu27'] == 1) & (selection(x))]) / len(df.loc[df['is_HLT_IsoMu27'] == 1])

    return np.vectorize(_diff), np.vectorize(_acc)

def diffacc2d(selection_x, selection_y):
    """
    selection_{x,y}:
        lumbda function that returns a bool like lambda x: df['pt'] > x
    return:
        function for difficiency depending on selections as defined in lambda function
        diff(x, y)
    """
    _diff = lambda x,y: len(df.loc[(df['is_HLT_IsoMu27'] == 1) & (selection_x(x)) & (selection_y(y))]) / len(df.loc[selection_x(x) & (selection_y(y))])
    _acc = lambda x,y: len(df.loc[(df['is_HLT_IsoMu27'] == 1) & (selection_x(x)) & (selection_y(y))]) / len(df.loc[df['is_HLT_IsoMu27'] == 1])

    return np.vectorize(_diff), np.vectorize(_acc)

diff_pt, acc_pt = diffacc(lambda x: df['pt'] > x)
#diff_eta, acc_eta = diffacc(lambda x: df['eta'] < x)

#diff_pfIso, acc_pfIso = diffacc(lambda x: df['pfIso'] < x)
#diff_tkIso, acc_tkIso = diffacc(lambda x: df['tkIso'] < x)

diff_tkIso_pfIso, acc_tkIso_pfIso = diffacc2d(lambda x: df['tkIso'] < x, lambda x: df['pfIso'] < x)
#diff_dxy_dz, acc_dxy_dz = diffacc2d(lambda x: df['dxy'] < x, lambda x: df['dz'] < x)

############
# plotting #
############

print(">>> plot")

def plot(range, diff, acc, xlabel='', plot_max_diff=False):
    xx = np.linspace(range[0], range[1], 200)
    plt.cla()
    yy_diff = diff(xx)
    yy_acc = acc(xx)

    plt.plot(xx, yy_diff, label='$\delta X$')
    #plt.plot(xx, yy_acc, label='$A$')

    if plot_max_diff:
        diff_max = xx[np.argmax(yy_diff)]
        ymin = min(yy_diff)#min(min(yy_diff), min(yy_acc))
        ymax = max(yy_diff)#max(max(yy_diff), max(yy_acc))
        plt.plot([diff_max,diff_max], [ymin,ymax], '-k', label='max($\delta X$)')

    plt.xlabel(xlabel)
    plt.legend()
    plt.savefig(output+"/triggerdiff_{0}_{1}to{2}.png".format(xlabel, str(range[0]).replace(".",""), str(range[1]).replace(".","")))

def plot2d(xrange, yrange, diff, acc, xlabel='', ylabel='', suffix=''):

    def _plot2d(func, zlabel, name):
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
            name, suffix))

    _plot2d(diff, '$\delta X$', 'difference')
    _plot2d(acc, '$A$', 'acceptance')

pdb.set_trace()
df = df.query('tkIso < 0.05 & pfIso < 0.12')

plot((27, 35), diff_pt, acc_pt, xlabel="pt", plot_max_diff=True)
#plot((27, 30), diff_pt, acc_pt, xlabel="pt")

#plot((27, 33), diff_pt, acc_pt, xlabel="pt")
#
#
# # plot((0.5, 2.4), diff_pt, acc_pt, xlabel="eta")
# # plot((0.025, 0.2), diff_tkIso, acc_tkIso, xlabel="tkIso")
# # plot((0.06, 0.4), diff_pfIso, acc_pfIso, xlabel="pfIso")
# #
# # plot((0.01, 1), diff_tkIso, acc_tkIso, xlabel="tkIso")
# # plot((0.01, 1), diff_pfIso, acc_pfIso, xlabel="pfIso")
#
# plot2d((0.025, 0.2), (0.06, 0.4), diff_tkIso_pfIso, acc_tkIso_pfIso, 'tkIso', 'pfIso')
#
# df = df.query('pt > 28')
# plot2d((0.025, 0.2), (0.06, 0.4), diff_tkIso_pfIso, acc_tkIso_pfIso, 'tkIso', 'pfIso','_pt28')
#
# df = df.query('pt > 30')
# plot2d((0.025, 0.2), (0.06, 0.4), diff_tkIso_pfIso, acc_tkIso_pfIso, 'tkIso', 'pfIso','_pt30')
#
# df = df.query('pt > 32')
# plot2d((0.025, 0.2), (0.06, 0.4), diff_tkIso_pfIso, acc_tkIso_pfIso, 'tkIso', 'pfIso','_pt32')
#
# #plot2d((0.1, 1), (0.25, 2), diff_dxy_dz, acc_dxy_dz, 'dxy', 'dz')
