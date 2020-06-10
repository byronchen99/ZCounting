from __future__ import division, print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from root_numpy import root2array, list_trees
from Utils import plot_scatter
from scipy.stats import pearsonr
import pdb

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df

import argparse

# plt.rcParams.update({'font.size': 18})

parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-i', '--input', nargs='+',
    help='specify input root file'
)
parser.add_argument(
    '-o', '--output', default='./',
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
branches = ['nPU', 'nPV',
            'muon_dxy', 'antiMuon_dxy',
            'muon_dz', 'antiMuon_dz',
            'muon_genVtxToPV', 'antiMuon_genVtxToPV',
            'muon_isFromPV', 'antiMuon_isFromPV'
            ]

print(">>> Load Events in gen acceptance")
df = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
df = pd.concat(df)

### dz vs PU

def plotMuonOverPU(quantity, xmax):
    dLowPU = np.concatenate(df.query("nPU<20")[["muon_{0}".format(quantity),"antiMuon_{0}".format(quantity)]].values)
    dMediumPU = np.concatenate(df.query("nPU>20 & nPU<45")[["muon_{0}".format(quantity),"antiMuon_{0}".format(quantity)]].values)
    dHighPU = np.concatenate(df.query("nPU>45")[["muon_{0}".format(quantity),"antiMuon_{0}".format(quantity)]].values)

    plt.clf()
    fig, ax = plt.subplots()

    # include overflow events bins in histogram
    dLowPU = map(lambda x: x if x<xmax else xmax-0.00001, dLowPU)
    dMediumPU =  map(lambda x: x if x<xmax else xmax-0.00001, dMediumPU)
    dHighPU =  map(lambda x: x if x<xmax else xmax-0.00001, dHighPU)

    ax.hist(dLowPU, bins=np.linspace(0, xmax, 20), label='nPU < 20', histtype='step', normed=True)
    ax.hist(dMediumPU, bins=np.linspace(0, xmax, 20), label='20 < nPU < 45', histtype='step', normed=True)
    ax.hist(dHighPU, bins=np.linspace(0, xmax, 20), label='nPU > 45', histtype='step', normed=True)

    ax.set(
        xlabel=quantity+" in cm",
        ylabel="normalized events",
        yscale='log',
    #    ylim=(min(hist-histErr)*0.99, max(hist+histErr)*1.01),
        xlim=(0., xmax)
           )
    ax.legend(loc=1)
    textstr = '\n'.join((
            r'$\mathrm{\mathbf{CMS}}$ Simulation $Z \rightarrow\ \mu\mu$',
            r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$',
            r'$p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4$',
        ))
    #box = dict(boxstyle='round', linecolor='white', facecolor='white', alpha=0.9)
    ax.text(0.1, 0.8, textstr, verticalalignment='bottom', transform=ax.transAxes) #, bbox=box

    plt.savefig(output+"/evntsVs{0}.png".format(quantity))
    plt.savefig(output+"/evntsVs{0}.pdf".format(quantity))

plotMuonOverPU("dz", 1.)
plotMuonOverPU("dxy", 0.4)

### avgDz over PU
xmin, xmax, xstep = 0, 70, 5
bin_edges = np.arange(xmin, xmax+xstep, xstep)
bin_centers = bin_edges[1:] - (bin_edges[1:] - bin_edges[:-1])/2
avgDz = []
avgDzUp = []
avgDzDown = []
avgDxy = []
avgDxyUp = []
avgDxyDown = []
for i in xrange(len(bin_edges) - 1):
    col = np.concatenate(df.query("nPU > {0} & nPU < {1}".format(bin_edges[i], bin_edges[i+1]))[["muon_dz","antiMuon_dz"]].values)
    avgDz.append(np.median(col))
    avgDzUp.append(np.percentile(col, 45))
    avgDzDown.append(np.percentile(col, 55))

    col = np.concatenate(df.query("nPU > {0} & nPU < {1}".format(bin_edges[i], bin_edges[i+1]))[["muon_dxy","antiMuon_dxy"]].values)
    avgDxy.append(np.median(col))
    avgDxyUp.append(np.percentile(col, 45))
    avgDxyDown.append(np.percentile(col, 55))


def plotAvgD(avgCol, upCol, downCol, name):
    avgCol = np.array(avgCol)
    upCol = np.array(upCol)
    downCol = np.array(downCol)

    plt.clf()
    fig, ax = plt.subplots()

    ax.fill_between(bin_centers, upCol, downCol, color='black', alpha=0.5, label='45% - 55% quantiles')
    ax.plot(bin_centers, avgCol, 'ko-', color='black', label='median')

    ran = max(avgCol) - min(avgCol)
    ax.set(
        xlabel="nPU",
        ylabel=name,
        #yscale='log'
        ylim=(min(avgCol)-ran*0.05, max(avgCol)+ran*0.05)
        )
    ax.legend(loc=4)

    textstr = '\n'.join((
            r'$\mathrm{\mathbf{CMS}}$ Simulation $Z \rightarrow\ \mu\mu$',
            r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$',
            r'$p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4$',
        ))
    ax.text(0.05, 0.8, textstr, verticalalignment='bottom', transform=ax.transAxes)

    plt.savefig(output+"/{0}VsPu.png".format(name))
    plt.savefig(output+"/{0}VsPu.pdf".format(name))

plotAvgD(avgDz, avgDzUp, avgDzDown, "dz")
plotAvgD(avgDxy,  avgDxyUp, avgDxyDown, "dxy")

exit()


### Vertex
evtFromPV = df.query('muon_isFromPV == 1 & antiMuon_isFromPV == 1')['nPV'].values
evtNotFromPV = df['nPV'].values

xmin, xmax, xstep = 0, 70, 5
bin_edges = np.arange(xmin, xmax+xstep, xstep)
bin_centers = bin_edges[1:] - (bin_edges[1:] - bin_edges[:-1])/2
bin_widths = bin_edges[1:] - bin_edges[:-1]
hist1, _ = np.histogram(evtFromPV, bins=bin_edges)
hist2, _ = np.histogram(evtNotFromPV, bins=bin_edges)

hist = hist1.astype('float64')/hist2
histErr = hist*np.sqrt(1./hist1 + 1./hist2)
fig, ax = plt.subplots()

ax.errorbar(bin_centers, hist, yerr=histErr, color='black')

ax.set(
    xlabel="nPV",
    ylabel="Events(PV) / Events",
    ylim=(min(hist-histErr)*0.99, max(hist+histErr)*1.01),
    xlim=(xmin, xmax)
       )

textstr = '\n'.join((
        r'$\mathrm{\mathbf{CMS}}$ Simulation $Z \rightarrow\ \mu\mu$',
        r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$',
        r'$p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4$',
    ))
box = dict(boxstyle='round', facecolor='white', alpha=0.9)
ax.text(0.05, 0.05, textstr, verticalalignment='bottom', bbox=box, transform=ax.transAxes)

plt.savefig(output+"/evntsVsNPV.png")
plt.savefig(output+"/evntsVsNPV.pdf")
