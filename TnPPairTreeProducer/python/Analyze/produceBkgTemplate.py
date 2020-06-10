from __future__ import division, print_function

# production of MC templates from TnP trees
#   Kernel Density Estimation (KDE) is used to extract the pdf shape
#   The pdf is then translated in a histogram which can be used for the fitting

import numpy as np
import pandas as pd
import os
import root_numpy as rn
import pdb
from ROOT import TH1D
import ROOT
import argparse
import matplotlib.pyplot as plt
import glob

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df

parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output dir'
)


args = parser.parse_args()

output = args.output


#ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/calculateDataEfficiency.C")
# ROOT.gROOT.SetBatch(True)

if os.path.isdir(output):
    print("output dir already exists, please remove or specify another name")
    exit()
os.mkdir(output)

MassMin_ = 56.
MassMax_ = 116.
MassBin_ = int(MassMax_ - MassMin_)

ptCut = 30.
tkIsoCut = 99999
dxyCut = 99999 #0.2
dzCut = 99999 #0.5

# acceptance selection
selection = 'pt1 >= {0} & pt2 >= {0} & q1 == q2'.format(ptCut)

if tkIsoCut is not None:
    selection += ' & (tkIso1 < {0} | (is2HLT & tkIso2 < {0}))'.format(tkIsoCut)
if dxyCut is not None:
    selection += ' & (dxy1 < {0} | (is2HLT & dxy2 < {0}))'.format(dxyCut)
if dzCut is not None:
    selection += ' & (dz1 < {0} | (is2HLT & dz2 < {0}))'.format(dzCut)

# specify which branches to load
branches = ['dilepMass', 'nPV', 'run', 'ls',
            'eta1', 'eta2',
            'pt1', 'pt2',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            'tkIso1', 'tkIso2',
            'dxy1', 'dxy2', 'dz1', 'dz2'
            ]

storage="/afs/desy.de/user/d/dwalter/store/SingleMuon/TnPPairTrees_V05/"

inputs = {
    ##'B':glob.glob(storage+"200324_212241/0000/output_0_*"),
    'C': glob.glob(storage+"200414_095025/0000/output_0_*"),
    'D': glob.glob(storage+"200414_095034/0000/output_0_*"),
    'E': glob.glob(storage+"200414_095047/0000/output_0_*"),
    ## 'F': glob.glob(storage+"200324_212331/0000/output_0_*"),
    'H':glob.glob(storage+"200415_062259/0000/output_0_*")
}

_df = []
for era, input in inputs.iteritems():
    print("##### era {0}".format(era))

    treeName = rn.list_trees(input[0])[0]
    print(">>> Load Events from {0} files".format(len(input)))
    _df.append(pd.concat([tree_to_df(rn.root2array(i, treeName, selection=selection, branches=branches)) for i in input]))


print(">>> Concatenate")
df = pd.concat(_df)
selected = df['is2HLT'] & (((df['tkIso1'] >= tkIsoCut) & (df['tkIso2'] < tkIsoCut))
                            | ((df['dxy1'] >= dxyCut) & (df['dxy2'] < dxyCut))
                            | ((df['dz1'] >= dzCut) & (df['dz2'] < dzCut))
                            )

to_switch = df[selected]
df = df[selected==False]
to_switch = to_switch.rename(index=str, columns={
    'tkIso1': 'tkIso2', 'tkIso2': 'tkIso1',
    'pt1': 'pt2', 'pt2': 'pt1',
    'eta1': 'eta2', 'eta2': 'eta1',
    'dxy1': 'dxy2', 'dxy2': 'dxy1',
    'dz1': 'dz2', 'dz2': 'dz1'
    })
df = pd.concat([df,to_switch],sort=True)
# downgrade is2HLT and isSel to isGlo
df['isGlo'] = df['isGlo'] \
    + df['is2HLT'] * (df['tkIso2'] >= tkIsoCut) + df['isSel'] * (df['tkIso2'] >= tkIsoCut) \
    + df['is2HLT'] * (df['dxy2'] >= dxyCut) + df['isSel'] * (df['dxy2'] >= dxyCut) \
    + df['is2HLT'] * (df['dz2'] >= dzCut) + df['isSel'] * (df['dz2'] >= dzCut)
df['is2HLT'] = df['is2HLT'] \
    * (df['tkIso2'] < tkIsoCut) \
    * (df['dxy2'] < dxyCut) \
    * (df['dz2'] < dzCut)
df['isSel'] = df['isSel'] \
    * (df['tkIso2'] < tkIsoCut) \
    * (df['dxy2'] < dxyCut) \
    * (df['dz2'] < dzCut)




from sklearn.neighbors.kde import KernelDensity

def analyzeShape(data, region, query=None):
    print('=== analyze shape for region {0}'.format(region))
    data = data.query(query) if query is not None else data
    xx = data['dilepMass'].values.reshape(-1, 1)
    xx_inBounds = data.query('dilepMass > {0} &  dilepMass < {1}'.format(MassMin_, MassMax_))['dilepMass'].values.reshape(-1, 1)

    MassDiff = MassMax_-MassMin_
    xgrid = np.linspace(MassMin_, MassMax_, 1000, endpoint=False).reshape(-1, 1)
#    xgrid_underflow = np.linspace(MassMin_-MassDiff, MassMin_, 1000, endpoint=False).reshape(-1, 1)
#    xgrid_overflow = np.linspace(MassMax_, MassMax_+MassDiff, 1000, endpoint=False).reshape(-1, 1)

    xgrid_hist = np.linspace(MassMin_+0.5, MassMax_-0.5, MassBin_).reshape(-1, 1)
#    xgrid_uf = np.linspace(MassMin_-MassDiff+0.5, MassMin_-0.5, MassBin_).reshape(-1, 1)
#    xgrid_of = np.linspace(MassMax_+0.5, MassMax_+MassDiff-0.5, MassBin_).reshape(-1, 1)

    results = []
    #kernels = ('gaussian', 'tophat', 'epanechnikov', 'exponential', 'linear', 'cosine')
    for krnl in ('gaussian', 'exponential', 'cosine'):
        print('fit {0} kernel'.format(krnl))

        min_nll = None
        bandwidths = np.logspace(0,2,10)#,base=2)
        #bandwidths = np.linspace(1, 10, 10)

        best_bw = bandwidths[0]

        for ibw, bw in enumerate(bandwidths):
            kde = KernelDensity(kernel=krnl, bandwidth=bw)

            def fit_kde(_xfit, _xval):
                _kde = kde.fit(_xfit)
                _xval = _xval[_xval>MassMin_]
                _xval = _xval[_xval<MassMax_]
                # negative log likelihood
                _nll = sum(-1 * _kde.score_samples(_xval.reshape(-1, 1)))
                return _nll

            nll = 0
            nfolds = 5

            xparts = np.array_split(xx[np.random.permutation(xx.shape[0])], nfolds)
            for k in range(nfolds): #k-fold cross validation fits
                xval = xparts[k]
                xfit = np.concatenate(xparts[:k] + xparts[k+1:])
                inll = fit_kde(xfit, xval)
                nll += inll

            if min_nll == None or min_nll > nll:
                min_nll = nll
                best_bw = bw

        print('best fit bw = {1} -> nll = {0}'.format(min_nll, best_bw))
        kde = KernelDensity(kernel=krnl, bandwidth=best_bw).fit(xx)
        ygrid = np.exp(kde.score_samples(xgrid))
        ygrid = ygrid/(sum(ygrid)*MassDiff/1000) * len(xx_inBounds)
         #() + np.exp(kde.score_samples(xgrid_underflow)[::-1]) + np.exp(kde.score_samples(xgrid_overflow)[::-1]))
        yy_hist = np.exp(kde.score_samples(xgrid_hist))
        yy_hist = yy_hist/(sum(yy_hist)*MassDiff/MassBin_) * len(xx_inBounds) #( + np.exp(kde.score_samples(xgrid_uf)[::-1]) + np.exp(kde.score_samples(xgrid_of)[::-1]))*len(xx)

        results.append([krnl, min_nll, best_bw, ygrid, yy_hist])

    idx_min = np.argmin([r[1] for r in results])

    # --- store as histogram
    print("store best histogram")
    tfile = ROOT.TFile.Open(output + "/{}.root".format(region), "RECREATE")

    for i, (krnl, mse, bw, ygrid, _) in enumerate(results):
        hBkgTemplate = TH1D("KDE_{0}".format(krnl), "th1d tnp same sign muons", MassBin_, MassMin_, MassMax_)
        #
        for i, y in enumerate(results[idx_min][-1]):
            hBkgTemplate.SetBinContent(i + 1, y)
        #
        hBkgTemplate.Write()
        hBkgTemplate.Delete()

    tfile.Close()

    # --- plotting
    print("plot")

    fig, ax = plt.subplots()
    plt.clf()
    #plt.title(krnl + ' kernal')
    plt.hist(data['dilepMass'], bins=np.linspace(MassMin_,MassMax_,MassBin_+1), histtype='stepfilled', color='0.5')
    for i, (krnl, mse, bw, ygrid, _) in enumerate(results):
        plt.plot(xgrid, ygrid, linewidth=3,
            label='{0}'.format(krnl+'*' if i == idx_min else krnl),
            #label='{0} (bw={1}, mse={2})'.format(krnl, bw,round(mse,2)),
            alpha=0.7)

    if 'central' in region:
        etaCut = '|\eta(\mu)| < 0.9'
    elif 'forward' in region:
        etaCut = '0.9 < |\eta(\mu)| < 2.4'
    else:
        etaCut = '|\eta(\mu)| < 2.4'

    box = dict(boxstyle='round', facecolor='white', alpha=0.9)
    textstr = '\n'.join((
        r'$\mathrm{{{0}}}$'.format(region.split("_")[0]),
        r'$p_\mathrm{{t}}(\mu) > {0}\ \mathrm{{GeV}} \qquad {1}$'.format(int(ptCut),etaCut),
        r'${0}\ \mathrm{{GeV}} < \mathrm{{M}}_{{\mu\mu}} < {1}\ \mathrm{{GeV}}$'.format(MassMin_,MassMax_),
        r'$\mathrm{{\mu}} = {0} \quad \mathrm{{\sigma}} = {1}$'.format(round(xx.mean(),1), round(xx.std(),1)),
    ))

    plt.text(0.02, 0.02, textstr, va='bottom', ha='left', bbox=box, linespacing=1.5, transform=ax.transAxes)
    plt.xlabel("Tag and probe mass")
    plt.ylabel("events")
    plt.xlim(MassMin_, MassMax_)
    #plt.ylim(len(data) / 70 * 0.5, len(data) / 70 * 1.5)
    plt.legend(loc=4)
    plt.savefig(output + '/kde_{0}.png'.format(region))


analyzeShape(df, 'Reco', 'is2HLT==1 | isSel==1')

analyzeShape(df, 'HLTPass_central', 'eta2 < 0.9 & is2HLT==1')
analyzeShape(df, 'HLTPass_forward', 'eta2 > 0.9 & is2HLT==1')
analyzeShape(df, 'HLTFail_central', 'eta2 < 0.9 & isSel==1')
analyzeShape(df, 'HLTFail_forward', 'eta2 > 0.9 & isSel==1')
analyzeShape(df, 'HLTPass', 'is2HLT==1')
analyzeShape(df, 'HLTFail', 'isSel==1')
#
analyzeShape(df, 'SelPass_central', 'eta2 < 0.9 & (is2HLT==1 | isSel==1)')
analyzeShape(df, 'SelPass_forward', 'eta2 > 0.9 & (is2HLT==1 | isSel==1)')
analyzeShape(df, 'SelFail_central', 'eta2 < 0.9 & isGlo==1')
analyzeShape(df, 'SelFail_forward', 'eta2 > 0.9 & isGlo==1')
analyzeShape(df, 'SelPass', '(is2HLT==1 | isSel==1)')
analyzeShape(df, 'SelFail', 'isGlo==1')
#
analyzeShape(df, 'GloPass_central', 'eta2 < 0.9 & (is2HLT==1 | isSel==1 | isGlo==1)')
analyzeShape(df, 'GloPass_forward', 'eta2 > 0.9 & (is2HLT==1 | isSel==1 | isGlo==1)')
analyzeShape(df, 'GloFail_central', 'eta2 < 0.9 & (isSta==1 | isTrk==1)')
analyzeShape(df, 'GloFail_forward', 'eta2 > 0.9 & (isSta==1 | isTrk==1)')
analyzeShape(df, 'GloPass_forward', '(is2HLT==1 | isSel==1 | isGlo==1)')
analyzeShape(df, 'GloFail_central', '(isSta==1 | isTrk==1)')
#
# analyzeShape(df, 'inclusive')
# analyzeShape(df, 'central', 'eta2 < 0.9')
# analyzeShape(df, 'forward', 'eta2 > 0.9')


# avgNPV = df['nPV'].mean()
# analyzeShape(df, 'nPVLow', 'nPV < {0}'.format(avgNPV))
# analyzeShape(df, 'nPVHigh', 'nPV > {0}'.format(avgNPV))
