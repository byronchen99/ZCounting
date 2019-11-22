from __future__ import division, print_function

# production of MC templates from TnP trees
#   Kernel Density Estimation (KDE) is used to extract the pdf shape
#   The pdf is then translated in a histogram which can be used for the fitting

import numpy as np
import pandas as pd
import os
from root_numpy import root2array, list_trees, array2hist
import pdb
from Utils.Utils import tree_to_df
from ROOT import TH1D
import ROOT
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-i', '--input', nargs='+',
    help='specify input tnp root files'
)
parser.add_argument(
    '--valData', nargs='+', default='',
    help='validation tnpFiles'
)
parser.add_argument(
    '--ptCut', type=float,
    help='specify lower pt cut on tag and probe muons'
)
parser.add_argument(
    '--kernel', type=str, default='gaussian',
    help='specify kernel function for KDE'
)
parser.add_argument(
    '--bw', type=float, default=5.,
    help='specify kernel bandwidth for KDE'
)
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output dir'
)


args = parser.parse_args()

inputs = args.input
output = args.output
valData = args.valData
ptCut = args.ptCut
kernel = args.kernel
bandwidth = args.bw

ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/calculateDataEfficiency.C")
# ROOT.gROOT.SetBatch(True)

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

# acceptance selection
selection = 'pt1 > {0} & pt2 > {0} & tkIso1/pt1<0.05'.format(ptCut)  # 'dilepMass > 66 ' \
# '& dilepMass < 116 ' \

# specify which branches to load
branches = ['dilepMass', 'nPV', # 'eventNumber', 'run', 'ls',
            'eta1', 'eta2',
            'pt1', 'pt2',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            'tkIso2',
            ]

print(">>> Load Events")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
print(">>> Concatenate")
df = pd.concat(df)

df['is2HLT'] = df['is2HLT'] * (df['tkIso2']/df['pt2'] < 0.05)
df['isSel'] = df['isSel'] + df['is2HLT'] * (df['tkIso2']/df['pt2'] > 0.05)

if valData != '':
    print(">>> Load Events")
    dfVal = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in valData]
    print(">>> Concatenate")
    dfVal = pd.concat(dfVal)
    dfVal['is2HLT'] = dfVal['is2HLT'] * (dfVal['tkIso2'] / dfVal['pt2'] < 0.05)
    dfVal['isSel'] = dfVal['isSel'] + dfVal['is2HLT'] * (dfVal['tkIso2'] / dfVal['pt2'] > 0.05)
else:
    dfVal = None

from sklearn.neighbors.kde import KernelDensity


def analyzeShape(data, region, query=None, vData=None):
    data = data.query(query) if query is not None else data
    xx = data['dilepMass'].values.reshape(-1, 1)
    xgrid = np.linspace(66., 116., 1000).reshape(-1, 1)

    for krnl in (kernel,): #('gaussian',):  # 'tophat', 'exponential', 'linear'):
        bw = [bandwidth, ]
        plt.clf()
        plt.title(krnl + ' kernal')
        plt.hist(data['dilepMass'], bins=70, histtype='stepfilled', color='0.5')
        yy = []
        for ibw in range(0, 1):
            print('fit {0} kernel with bw = {1}'.format(krnl, bw[ibw]))
            kde = KernelDensity(kernel=krnl, bandwidth=bw[ibw]).fit(xx)

            yy.append(np.exp(kde.score_samples(xgrid)))
            plt.plot(xgrid, len(data) * yy[ibw], linewidth=3, label='bw = {0}'.format(bw[ibw]),
                     alpha=0.7)
            bw.append(bw[ibw] + 2)
        plt.text(70., len(data) / 70 * 1.4, region)
        plt.xlabel("Tag and probe mass")
        plt.ylabel("events")
        plt.xlim(66., 116.)
        #plt.ylim(len(data) / 70 * 0.5, len(data) / 70 * 1.5)
        plt.legend(loc=1)
        plt.savefig(output + '/kde_{0}_{1}.png'.format(krnl, region))

        # validation plot
        if vData is not None:
            vData = vData.query(query) if query is not None else vData
            plt.clf()
            plt.title(krnl + ' kernal - validation')
            plt.hist(vData['dilepMass'], bins=50, histtype='stepfilled', color='0.5')
            for iy, y in enumerate(yy):
                plt.plot(xgrid, len(vData) * 1. / (50/1000. * sum(y)) * y, linewidth=3, label='bw = {0}'.format(bw[iy]),
                         alpha=0.7)
            plt.text(70., len(vData) / 50 * 1.4, region)
            plt.xlabel("Tag and probe mass")
            plt.ylabel("events")
            plt.xlim(66., 116.)
            plt.ylim(len(vData) / 50 * 0.5, len(vData) / 50 * 1.5)
            plt.legend(loc=1)
            plt.savefig(output + '/kde_val_{0}_{1}.png'.format(krnl, region))


analyzeShape(df, 'inclusive', vData=dfVal)
analyzeShape(df, 'central', 'eta2 < 0.9', dfVal)
analyzeShape(df, 'forward', 'eta2 > 0.9', dfVal)

analyzeShape(df, 'GloPass', 'is2HLT==1 | isSel==1 | isGlo==1', dfVal)
analyzeShape(df, 'GloFail', 'isSta==1 | isTrk==1', dfVal)
analyzeShape(df, 'SelPass', 'is2HLT==1 | isSel==1', dfVal)
analyzeShape(df, 'SelFail', 'isGlo==1', dfVal)
analyzeShape(df, 'HLTPass', 'is2HLT==1 ', dfVal)
analyzeShape(df, 'HLTFail', 'isSel==1', dfVal)

avgNPV = df['nPV'].mean()
analyzeShape(df, 'nPVLow', 'nPV < {0}'.format(avgNPV), dfVal)
analyzeShape(df, 'nPVHigh', 'nPV > {0}'.format(avgNPV), dfVal)

tfile = ROOT.TFile.Open(output + "/template_Pt{0}K{1}BW{2}.root".format(ptCut, kernel, bandwidth), "RECREATE")

def storeShape(data, region, query=None):
    print(">>> Create histogram shape in root file")
    data = data.query(query) if query is not None else data
    xx = data['dilepMass'].values.reshape(-1, 1)
    hBkgTemplate = TH1D("bkg_template_{0}".format(region), "th1d tnp same sign muons inclusive ", MassBin, MassMin,
                        MassMax)

    kde = KernelDensity(kernel=kernel, bandwidth=bandwidth).fit(xx)
    xgrid = np.linspace(MassMin+0.5, MassMax-0.5, MassBin).reshape(-1, 1)
    yy = np.exp(kde.score_samples(xgrid))

    for i, y in enumerate(yy):
        hBkgTemplate.SetBinContent(i + 1, y)

    hBkgTemplate.Write()
    hBkgTemplate.Delete()


storeShape(df, 'central', 'eta2 < 0.9')
storeShape(df, 'forward', 'eta2 > 0.9')
storeShape(df, 'inclusive', 'eta2 > 0.9')
storeShape(df, 'GloPass', 'is2HLT==1 | isSel==1 | isGlo==1')
storeShape(df, 'GloFail', 'isSta==1 | isTrk==1')
storeShape(df, 'SelPass', 'is2HLT==1 | isSel==1')
storeShape(df, 'SelFail', 'isGlo==1')
storeShape(df, 'HLTPass', 'is2HLT==1 ')
storeShape(df, 'HLTFail', 'isSel==1')

tfile.Close()
