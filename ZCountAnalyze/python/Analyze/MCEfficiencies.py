from __future__ import division, print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from root_numpy import root2array, list_trees
from Utils import tree_to_df
import math

import argparse

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


def eff(nReco, nAcc):
    """
    MC Z Effciciency
    :param nReco:
    :param nAcc:
    :param nBoth: intersect of nReco and nAcc for uncertainty estimation
    :return:
    """
    if nAcc < 10 or nReco < 5:
        return 0., 0.
    eff = nReco / nAcc
    x = nReco
    b = nAcc - nReco
    err = np.sqrt((b / (x + b) ** 2) ** 2 * x
                  + (x / (x + b) ** 2) ** 2 * b
                  )
    return eff, err


def zMCEff1D(df, bins, obs, region='inclusive', sel=''):
    """
    compute Z efficiency from MC and Z->mumu efficiency from tag and probe
    :param df: dataframe with events satisfying the gen accepance
    :param bins:
    :param obs:  oberservable
    :param region: string for naming reasons
    :param sel: selection on gen parameters
    :return:
    """
    print(">>> make differential efficiencies for " + obs + " in " + region + " region")

    if sel != '':
        df = df.query(sel)

    Efficiency.Set(bins, obs, region)

    for i in range(0, len(bins) - 1):
        bin_low = bins[i]
        bin_high = bins[i + 1]
        sel = '{0}>{1} & {0}<{2}'.format(obs, bin_low, bin_high)

        Efficiency.Calculate(df.query(sel), i)

    Efficiency.Plot()


class Efficiency:
    collection = []

    @classmethod
    def Set(cls, bins, obs, region):
        for eff in cls.collection:
            eff.set(bins, obs, region)

    @classmethod
    def Calculate(cls, df, ibin):
        for eff in cls.collection:
            eff.calculate(df, ibin)

    @classmethod
    def Plot(cls):
        for eff in cls.collection:
            eff.plot()

    def __init__(self, name, reqID, reqPFIso, reqTkIso, reqHLT=None, triggerName=''):
        self.name = name
        self.ID = reqID
        self.pfIso = reqPFIso
        self.tkIso = reqTkIso
        self.HLT = reqHLT
        self.triggerName = triggerName
        self.eff_true = []
        self.eff_tnp = []
        self.bins = []
        self.obs = None
        self.region = None
        self.x = []
        Efficiency.collection.append(self)

    def set(self, bins, obs, region):
        self.bins = bins
        self.obs = obs
        self.region = region
        self.x = bins[:-1] + (bins[1:] - bins[:-1]) / 2
        self.eff_true = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
        self.eff_tnp = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]

    def calculate(self, df, ibin):
        denom = len(df.query('muon_recoMatches == 1 & muon_ID>={0} & muon_pfIso >= {1} & muon_tkIso >= {2}'
                             .format(self.ID, self.pfIso, self.tkIso)))
        nom = len(df.query('muon_recoMatches == 1 & antiMuon_recoMatches == 1'
                           '& muon_ID >= {0} & antiMuon_ID >= {0}'
                           '& muon_pfIso >= {1} & antiMuon_pfIso >= {1}'
                           '& muon_tkIso >= {2} & antiMuon_tkIso >= {2}'
                           .format(self.ID, self.pfIso, self.tkIso)))
        MuEff, MuEff_Err = eff(nom, denom)

        if self.HLT:
            denom = len(df.query('muon_recoMatches == 1 & antiMuon_recoMatches == 1'
                                 '& muon_ID >= {0} & antiMuon_ID >= {0}'
                                 '& muon_pfIso >= {1} & antiMuon_pfIso >= {1}'
                                 '& muon_tkIso >= {2} & antiMuon_tkIso >= {2}'
                                 '& muon_hlt_{3} == 1'
                                 .format(self.ID, self.pfIso, self.tkIso, self.HLT)))
            nom = len(df.query('muon_ID>={0} & muon_recoMatches == 1'
                               '& antiMuon_ID>={0} & antiMuon_recoMatches == 1'
                               '& muon_pfIso >= {1} & muon_tkIso >= {2}'
                               '& antiMuon_pfIso >= {1} & antiMuon_tkIso >= {2}'
                               '& muon_hlt_{3} == 1 & antiMuon_hlt_{3} == 1'
                               .format(self.ID, self.pfIso, self.tkIso, self.HLT)))
            MuEff_HLT, MuEff_HLT_Err = eff(nom, denom)

            self.eff_tnp[0][ibin] = (1 - (1 - MuEff_HLT) ** 2) * MuEff ** 2
            self.eff_tnp[1][ibin] = np.sqrt((2 * (1 - MuEff_HLT) * MuEff ** 2 * MuEff_HLT_Err) ** 2 \
                                            + ((1 - (1 - MuEff_HLT) ** 2) * 2 * MuEff * MuEff_HLT_Err) ** 2)
        else:
            self.eff_tnp[0][ibin] = MuEff ** 2
            self.eff_tnp[1][ibin] = 2 * MuEff * MuEff_Err

        nZAcc = len(df)

        query = 'muon_recoMatches == 1 & antiMuon_recoMatches ==1 ' \
                '& muon_ID>={0} & antiMuon_ID>={0} ' \
                '& muon_pfIso >= {1} & antiMuon_pfIso >= {1} ' \
                '& muon_tkIso >= {2} & antiMuon_tkIso >= {2}'.format(self.ID, self.pfIso, self.tkIso)
        if self.HLT:
            query += '& ( muon_hlt_{0} == 1 | antiMuon_hlt_{0} == 1)'.format(self.HLT)
        nZReco = len(df.query(query))

        self.eff_true[0][ibin], self.eff_true[1][ibin] = eff(nZReco, nZAcc)

    def plot(self):
        plt.clf()
        fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        fig.subplots_adjust(hspace=0)

        ax[0].errorbar(self.x, self.eff_true[0], xerr=np.zeros(len(self.x)), yerr=self.eff_true[1], fmt='bo',
                       label='true')
        ax[0].errorbar(self.x, self.eff_tnp[0], xerr=np.zeros(len(self.x)), yerr=self.eff_tnp[1], fmt='ro', label='tnp')
        ymin = min([min(self.eff_true[0]), min(self.eff_tnp[0])])
        ymax = max([max(self.eff_true[0]), max(self.eff_tnp[0])])
        ymax += (ymax-ymin)*0.4
        ax[0].set(xlim=(self.bins[0], self.bins[-1]), ylim=(ymin, ymax))
        ax[0].legend()
        ax[0].text(0.05, 0.9, r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$',
                   transform=ax[0].transAxes)
        ax[0].text(0.05, 0.82, r'$p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4$',
                   transform=ax[0].transAxes)
        ax[0].text(0.05, 0.74, "${}$".format(self.region), transform=ax[0].transAxes)
        if self.HLT:
            ax[0].text(0.05, 0.04,
                       r'{0}'.format(self.triggerName),
                       transform=ax[0].transAxes, color='g')
        ax[0].set_ylabel(r'$\epsilon_\mathrm{Z}$')
        ax[0].set_xlabel(self.obs)
        # ax[0].set_yticks([0.8, 0.85, 0.9, 0.95, 1.0, 1.05])
        ax[0].set_title("Z efficiency at {0} muon level".format(self.name))

        pulls = self.eff_tnp[0] - self.eff_true[0]
        pulls_sig = np.sqrt(self.eff_tnp[1] ** 2 + self.eff_true[1] ** 2)

        ax[1].errorbar(self.x, pulls, xerr=np.zeros(len(self.x)), yerr=pulls_sig,
                       fmt='ko')  # , label='factorized - true')
        ax[1].plot(self.x, np.zeros(len(self.x)), color='gray', linestyle='dashed')
        ax[1].set_ylim(-0.01, 0.03)
        ax[1].set_ylabel(r'$\epsilon^\mathrm{tnp}_\mathrm{Z} - \epsilon^\mathrm{true}_\mathrm{Z}$')
        ax[1].set_xlabel(self.obs)
        ax[1].set_yticks([-0.01, 0., 0.01, 0.02])
        plt.savefig(output + '/ZMuMu_{0}_{1}_{2}_level.png'.format(self.obs, self.region, self.name))
        plt.close()


# acceptance selection
selection = 'z_genMass > 66 ' \
            '& z_genMass < 116 ' \
            '& muon_genPt > 27 ' \
            '& antiMuon_genPt > 27 ' \
            '& abs(muon_genEta) < 2.4 ' \
            '& abs(antiMuon_genEta) < 2.4 ' \
            '& muon_recoMatches <= 1' \
            '& antiMuon_recoMatches <= 1'

# specify which branches to load
branches = ['nPV', 'nPU', 'z_recoMass', 'z_genMass',
            'muon_genEta', 'antiMuon_genEta',
            'muon_genPhi', 'antiMuon_genPhi',
            'muon_genPt', 'antiMuon_genPt',
            'muon_recoMatches', 'antiMuon_recoMatches',
            'muon_ID', 'antiMuon_ID',
            'muon_triggerBits', 'antiMuon_triggerBits',
            'muon_tkIso', 'antiMuon_tkIso',
            'muon_pfIso', 'antiMuon_pfIso'
            ]

print(">>> Load Events in gen acceptance")
dfGen = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
dfGen = pd.concat(dfGen)

print(">>> add new columns")
dfGen['delPhiLL'] = abs(
    abs(dfGen['muon_genPhi'] - dfGen['antiMuon_genPhi']).apply(lambda x: x - 2 * math.pi if x > math.pi else x))
dfGen['delEtaLL'] = abs(dfGen['muon_genEta'] - dfGen['antiMuon_genEta'])
dfGen['delRLL'] = np.sqrt(
    (dfGen['muon_genEta'] - dfGen['antiMuon_genEta']) ** 2 + (dfGen['muon_genPhi'] - dfGen['antiMuon_genPhi']) ** 2)
dfGen['delPtLL'] = abs(dfGen['muon_genPt'] - dfGen['antiMuon_genPt'])
dfGen['relPtLL'] = abs(dfGen['muon_genPt'] - dfGen['antiMuon_genPt']) / abs(
    dfGen['muon_genPt'] + dfGen['antiMuon_genPt'])
dfGen['sumPtLL'] = dfGen['muon_genPt'] + dfGen['antiMuon_genPt']

dfGen = dfGen.query('delRLL > 0.4')

print(">>> convert bit code into bit map")
for iBit in range(0, 5):
    nBit = 2 ** iBit
    iBit += 1
    dfGen['muon_hlt_{0}'.format(iBit)] = dfGen['muon_triggerBits'].apply(
        lambda x: 1 if x % (nBit * 2) >= nBit else 0)
    dfGen['antiMuon_hlt_{0}'.format(iBit)] = dfGen['antiMuon_triggerBits'].apply(
        lambda x: 1 if x % (nBit * 2) >= nBit else 0)

# --- define Efficiencies
Efficiency("ZNoID", 0, 0, 0)
Efficiency("ZTrk", 1, 0, 0)
Efficiency("ZTrkOrSta", 2, 0, 0)
Efficiency("ZGlobal", 3, 0, 0)
Efficiency("ZcTight", 4, 0, 0)
Efficiency("ZTight", 5, 0, 0)
Efficiency("ZcTightL1SMu18", 4, 0, 0, 1, '$\mathrm{HLT\_L1SingleMu18\_v^{*}}$')
Efficiency("ZcTightL1SMu25", 4, 0, 0, 2, '$\mathrm{HLT\_L1SingleMu25\_v^{*}}$')
Efficiency("ZcTightIsoMu24", 4, 0, 0, 3, '$\mathrm{HLT\_IsoMu24\_v^{*}}$')
Efficiency("ZcTightIsoMu27", 4, 0, 0, 4, '$\mathrm{HLT\_IsoMu27\_v^{*}}$')
Efficiency("ZcTightIsoMu30", 4, 0, 0, 5, '$\mathrm{HLT\_IsoMu30\_v^{*}}$')

zMCEff1D(dfGen, np.linspace(0.5, 74.5, 25), 'nPU', 'inclusive')
zMCEff1D(dfGen, np.linspace(0.5, 74.5, 25), 'nPU', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPU', 'BE',
#       sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(0.5, 74.5, 25), 'nPU', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

# --- differential efficiencies
# zMCEff1D(dfGen, np.linspace(0., 5., 20), 'delRLL')
# zMCEff1D(dfGen, np.linspace(0., 5., 20), 'delRLL', 'BB',
#          sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
# zMCEff1D(dfGen, np.linspace(0., 5., 20), 'delRLL', 'BE',
#          sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
# zMCEff1D(dfGen, np.linspace(0., 5., 20), 'delRLL', 'EE',
#          sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')
#
# zMCEff1D(dfGen, np.linspace(0., 2.5, 20), 'delEtaLL')
# zMCEff1D(dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'BB',
#          sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
# zMCEff1D(dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'BE',
#          sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
# zMCEff1D(dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'EE',
#          sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')
#
# zMCEff1D(dfGen, np.linspace(0., math.pi, 20), 'delPhiLL')
# zMCEff1D(dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'BB',
#          sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
# zMCEff1D(dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'BE',
#          sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
# zMCEff1D(dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'EE',
#          sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')
#
# zMCEff1D(dfGen, np.linspace(0., 100, 20), 'delPtLL')
# zMCEff1D(dfGen, np.linspace(0., 100, 20), 'delPtLL', 'BB',
#          sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
# zMCEff1D(dfGen, np.linspace(0., 100, 20), 'delPtLL', 'BE',
#          sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
# zMCEff1D(dfGen, np.linspace(0., 100, 20), 'delPtLL', 'EE',
#          sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')
#
# zMCEff1D(dfGen, np.linspace(0., 1., 20), 'relPtLL')
# zMCEff1D(dfGen, np.linspace(0., 1., 20), 'relPtLL', 'BB',
#          sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
# zMCEff1D(dfGen, np.linspace(0., 1., 20), 'relPtLL', 'BE',
#          sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
# zMCEff1D(dfGen, np.linspace(0., 1., 20), 'relPtLL', 'EE',
#          sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')
#
# zMCEff1D(dfGen, np.linspace(50., 250, 20), 'sumPtLL')
# zMCEff1D(dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'BB',
#          sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
# zMCEff1D(dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'BE',
#          sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
# zMCEff1D(dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'EE',
#          sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')
#
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'inclusive')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'BB',
#          sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'BE',
#          sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'EE',
#          sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delRLL0to25', sel='delRLL < 2.5')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delRLL25toInf', sel='delRLL > 2.5')
#
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delEtaLL0to1', sel='delEtaLL < 1')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delEtaLL1toInf', sel='delEtaLL > 1')
#
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPhiLL0to25', sel='delPhiLL < 2.5')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPhiLL25toPi', sel='delPhiLL > 2.5')
#
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPtLL0to10', sel='delPtLL < 15')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPtLL10toInf', sel='delPtLL > 15')
#
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'relPtLL0to01', sel='relPtLL < 0.1')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'relPtLL01toInf', sel='relPtLL > 0.1')
#
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPU', 'sumPtLL0to100', sel='sumPtLL < 100')
# zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPU', 'sumPtLL100toInf', sel='sumPtLL > 100')
