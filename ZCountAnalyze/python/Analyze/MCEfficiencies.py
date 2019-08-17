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

    ZEff = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    ZEff_sel = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    ZEff_glo = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]

    MuEff_HLT = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Sel = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Glo = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]

    ZMuMuEff = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    ZMuMuEff_sel = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    ZMuMuEff_glo = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]

    for i in range(0, len(bins) - 1):
        bin_low = bins[i]
        bin_high = bins[i + 1]
        sel = '{0}>{1} & {0}<{2}'.format(obs, bin_low, bin_high)

        ZReco_i = df.query(sel)

        # total glo muon efficiency
        denomGlo = len(ZReco_i.query('muon_Category>=3 & muon_recoMatches == 1'))
        nomGlo = len(ZReco_i.query('muon_Category>=3 & muon_recoMatches == 1'
                                   '& antiMuon_Category>=3 & antiMuon_recoMatches ==1'))
        MuEff_Glo[0][i], MuEff_Glo[1][i] = eff(nomGlo, denomGlo)

        # total sel efficiency
        denomSel = len(ZReco_i.query('muon_Category>=4 & muon_recoMatches == 1'))
        nomSel = len(ZReco_i.query('muon_Category>=4 & muon_recoMatches == 1'
                                   '& antiMuon_Category>=4 & antiMuon_recoMatches ==1'))
        MuEff_Sel[0][i], MuEff_Sel[1][i] = eff(nomSel, denomSel)

        ZMuMuEff_glo[0][i] = MuEff_Glo[0][i] ** 2
        ZMuMuEff_glo[1][i] = 2 * MuEff_Glo[0][i] * MuEff_Glo[1][i]
        ZMuMuEff_sel[0][i] = MuEff_Sel[0][i] ** 2
        ZMuMuEff_sel[1][i] = 2 * MuEff_Sel[0][i] * MuEff_Sel[1][i]

        # hlt efficiency in respect of Sel muons
        denomHLT = len(ZReco_i.query('muon_Category>=5 & muon_recoMatches == 1'
                                     '& antiMuon_Category>=4 & antiMuon_recoMatches == 1'))
        nomHLT = len(ZReco_i.query('muon_Category>=5 & muon_recoMatches == 1'
                                   '& antiMuon_Category>=5 & antiMuon_recoMatches == 1'))
        MuEff_HLT[0][i], MuEff_HLT[1][i] = eff(nomHLT, denomHLT)

        ZMuMuEff[0][i] = (1 - (1 - MuEff_HLT[0][i]) ** 2) * MuEff_Sel[0][i] ** 2
        ZMuMuEff[1][i] = np.sqrt((2 * (1 - MuEff_HLT[0][i]) * MuEff_Sel[0][i] ** 2 * MuEff_HLT[1][i]) ** 2 \
                                 + ((1 - (1 - MuEff_HLT[0][i]) ** 2) * 2 * MuEff_Sel[0][i] * MuEff_Sel[1][i]) ** 2)

        #  full true Z efficiency
        nZAcc_i = len(df.query(sel))
        nZReco_i = len(df.query(sel + '& muon_recoMatches == 1 & antiMuon_recoMatches ==1'
                                      '& (muon_Category==5 | antiMuon_Category==5) & muon_Category>=4 & antiMuon_Category>=4'))
        ZEff[0][i], ZEff[1][i] = eff(nZReco_i, nZAcc_i)

        #  true Z efficiency at selection level (without requiring trigger)
        nZAcc_i = len(df.query(sel))
        nZReco_i = len(df.query(sel + '& muon_recoMatches == 1 & antiMuon_recoMatches ==1 '
                                      '& muon_Category>=4 & antiMuon_Category>=4'))
        ZEff_sel[0][i], ZEff_sel[1][i] = eff(nZReco_i, nZAcc_i)

        #  true Z efficiency at gobal muon level (requiring two global muons)
        nZAcc_i = len(df.query(sel))
        nZReco_i = len(df.query(sel + '& muon_recoMatches == 1 & antiMuon_recoMatches ==1 '
                                      '& muon_Category>=3 & antiMuon_Category>=3'))
        ZEff_glo[0][i], ZEff_glo[1][i] = eff(nZReco_i, nZAcc_i)

    x = bins[:-1] + (bins[1:] - bins[:-1]) / 2

    def plot_zeff(eff_true, eff_tnp, name):
        plt.clf()
        fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        fig.subplots_adjust(hspace=0)

        ax[0].errorbar(x, eff_true[0], xerr=np.zeros(len(x)), yerr=eff_true[1], fmt='bo', label='true')
        ax[0].errorbar(x, eff_tnp[0], xerr=np.zeros(len(x)), yerr=eff_tnp[1], fmt='ro', label='tnp')
        ax[0].set(xlim=(bins[0], bins[-1]), ylim=(0.75, 1.1))
        ax[0].legend()
        ax[0].text(0.05, 0.9, r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$',
                   transform=ax[0].transAxes)
        ax[0].text(0.05, 0.82, r'$p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4$',
                   transform=ax[0].transAxes)
        ax[0].text(0.05, 0.74, region, transform=ax[0].transAxes)
        ax[0].set_ylabel(r'$\epsilon_\mathrm{Z}$')
        ax[0].set_xlabel(obs)
        ax[0].set_yticks([0.8, 0.85, 0.9, 0.95, 1.0, 1.05])
        ax[0].set_title("Z efficiency at {0} muon level".format(name))

        pulls = eff_tnp[0] - eff_true[0]
        pulls_sig = np.sqrt(eff_tnp[1] ** 2 + eff_true[1] ** 2)

        ax[1].errorbar(x, pulls, xerr=np.zeros(len(x)), yerr=pulls_sig, fmt='ko')  # , label='factorized - true')
        ax[1].plot(x, np.zeros(len(x)), color='gray', linestyle='dashed')
        ax[1].set_ylim(-0.02, 0.04)
        ax[1].set_ylabel(r'$\epsilon^\mathrm{tnp}_\mathrm{Z} - \epsilon^\mathrm{true}_\mathrm{Z}$')
        ax[1].set_xlabel(obs)
        ax[1].set_yticks([-0.01, 0., 0.01, 0.02, 0.03])
        plt.savefig(output + '/ZMuMu_{0}_{1}_{2}_level.png'.format(obs, region, name))
        plt.close()

    plot_zeff(ZEff, ZMuMuEff, "mediumSelHLT")
    plot_zeff(ZEff_sel, ZMuMuEff_sel, "mediumSel")
    plot_zeff(ZEff_glo, ZMuMuEff_glo, "global")


# acceptance selection
selection = 'z_genMass > 66 ' \
            '& z_genMass < 116 ' \
            '& muon_genPt > 27 ' \
            '& antiMuon_genPt > 27 ' \
            '& abs(muon_genEta) < 2.4 ' \
            '& abs(antiMuon_genEta) < 2.4 '

# specify which branches to load
branches = ['nPV', 'nPU', 'z_recoMass', 'z_genMass',
            'muon_genEta', 'antiMuon_genEta',
            'muon_genPhi', 'antiMuon_genPhi',
            'muon_genPt', 'antiMuon_genPt',
            'muon_recoMatches', 'antiMuon_recoMatches',
            'muon_CategoryMedium', 'antiMuon_CategoryMedium',
            'muon_CategoryTight', 'antiMuon_CategoryTight'
            ]

print(">>> Load Events in gen acceptance")
dfGen = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches), 5) for i in inputs]
dfGen = pd.concat(dfGen)

dfGen = dfGen.rename(columns={"muon_CategoryMedium": "muon_Category", "antiMuon_CategoryMedium": "antiMuon_Category"})
#dfGen = dfGen.rename(columns={"muon_CategoryTight": "muon_Category", "antiMuon_CategoryTight": "antiMuon_Category"})

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

# --- differential efficiencies
dfList = []
zMCEff1D(dfGen, np.linspace(0., 5., 20), 'delRLL')
zMCEff1D(dfGen, np.linspace(0., 5., 20), 'delRLL', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
zMCEff1D(dfGen, np.linspace(0., 5., 20), 'delRLL', 'BE',
         sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(0., 5., 20), 'delRLL', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

zMCEff1D(dfGen, np.linspace(0., 2.5, 20), 'delEtaLL')
zMCEff1D(dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
zMCEff1D(dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'BE',
         sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

zMCEff1D(dfGen, np.linspace(0., math.pi, 20), 'delPhiLL')
zMCEff1D(dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
zMCEff1D(dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'BE',
         sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

zMCEff1D(dfGen, np.linspace(0., 100, 20), 'delPtLL')
zMCEff1D(dfGen, np.linspace(0., 100, 20), 'delPtLL', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
zMCEff1D(dfGen, np.linspace(0., 100, 20), 'delPtLL', 'BE',
         sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(0., 100, 20), 'delPtLL', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

zMCEff1D(dfGen, np.linspace(0., 1., 20), 'relPtLL')
zMCEff1D(dfGen, np.linspace(0., 1., 20), 'relPtLL', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
zMCEff1D(dfGen, np.linspace(0., 1., 20), 'relPtLL', 'BE',
         sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(0., 1., 20), 'relPtLL', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

zMCEff1D(dfGen, np.linspace(50., 250, 20), 'sumPtLL')
zMCEff1D(dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
zMCEff1D(dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'BE',
         sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'inclusive')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'BE',
         sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPU', 'inclusive')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPU', 'BB',
         sel='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPU', 'BE',
         sel='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPU', 'EE',
         sel='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9')

zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delRLL0to25', sel='delRLL < 2.5')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delRLL25toInf', sel='delRLL > 2.5')

zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delEtaLL0to1', sel='delEtaLL < 1')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delEtaLL1toInf', sel='delEtaLL > 1')

zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPhiLL0to25', sel='delPhiLL < 2.5')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPhiLL25toPi', sel='delPhiLL > 2.5')

zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPtLL0to10', sel='delPtLL < 15')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPtLL10toInf', sel='delPtLL > 15')

zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'relPtLL0to01', sel='relPtLL < 0.1')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'relPtLL01toInf', sel='relPtLL > 0.1')

zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'sumPtLL0to100', sel='sumPtLL < 100')
zMCEff1D(dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'sumPtLL100toInf', sel='sumPtLL > 100')
