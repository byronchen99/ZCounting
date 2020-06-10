from __future__ import division, print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from root_numpy import root2array, list_trees
import math

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df

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

if(len(treeName) > 1):
    print("more then one tree in file ... specify, which tree to use")
    exit()


def muTagAndProbeEff(hlt, pp, fp):
    """
    Muon tag and probe efficiency
    :param hlt: passing probes which also pass trigger selection
    :param pp:  other passing probes
    :param fp:  failing probes
    :return: eff, eff_err:
    """
    if 2*hlt+pp <= 0:
        return 0., 0.
    elif fp <= 0:
        return 0., 0.

    eff = (2*hlt+pp)/(2*hlt+pp+fp)
    eff_err = 1. / (2*hlt+fp+pp)**2 * np.sqrt(pp*fp**2 + hlt*(2*fp)**2 + fp*(pp+2*hlt)**2)
    return eff, eff_err


def eff_ztomumu_old(nHLT_,nSel_,nGlo_,nSta_,nTrk_):
    eff_HLT, err_HLT = muTagAndProbeEff(nHLT_, 0, nSel_)
    eff_Sel, err_Sel = muTagAndProbeEff(nHLT_, nSel_, nGlo_)
    eff_Glo, err_Glo = muTagAndProbeEff(nHLT_, nSel_ + nGlo_, nTrk_ + nSta_)

    if eff_HLT <= 0 or eff_Sel <=0 or eff_Glo <= 0:
        return 0., 0.

    eff_tot = (1 - (1 - eff_HLT)**2) * eff_Sel**2 * eff_Glo**2
    err_tot = 2 * eff_tot * np.sqrt(((1-eff_HLT)/(1-(1-eff_HLT)**2)*err_HLT)**2
                                    + (err_Sel/eff_Sel)**2 + (err_Glo/eff_Glo)**2
                                    )
    return eff_tot, err_tot


def eff_ztomumu_new(nHLT_, nSel_, nGlo_, nSta_, nTrk_):
    eff_HLT, err_HLT = muTagAndProbeEff(nHLT_, 0, nSel_)
    eff_Sel, err_Sel = muTagAndProbeEff(nHLT_, nSel_, nGlo_)
    eff_Trk, err_Trk = muTagAndProbeEff(nHLT_, nSel_ + nGlo_ + nTrk_, nSta_)
    eff_Sta, err_Sta = muTagAndProbeEff(nHLT_, nSel_ + nGlo_ + nSta_, nTrk_)

    if eff_HLT <= 0 or eff_Sel <= 0 or eff_Trk <= 0 or eff_Sta <= 0:
        return 0., 0.

    eff_tot = (1 - (1 - eff_HLT)**2) * eff_Sel**2 * eff_Trk**2 * eff_Sta**2
    err_tot = 2 * eff_tot * np.sqrt(((1-eff_HLT)/(1-(1-eff_HLT)**2)*err_HLT)**2
                                    + (err_Sel/eff_Sel)**2 + (err_Trk/eff_Trk)**2 + (err_Sta/eff_Sta)**2
                                    )
    return eff_tot, err_tot


def eff_ztomumu_sel(nHLT_, nSel_, nGlo_, nSta_, nTrk_):
    #  facctorised Z efficiency from tag and probe at selection level
    eff_Sel, err_Sel = muTagAndProbeEff(nHLT_, nSel_, nGlo_)
    eff_Glo, err_Glo = muTagAndProbeEff(nHLT_, nSel_ + nGlo_, nTrk_ + nSta_)

    if eff_Sel <= 0 or eff_Glo <= 0:
        return 0., 0.

    eff_tot = eff_Sel ** 2 * eff_Glo ** 2
    err_tot = 2 * eff_tot * np.sqrt((err_Sel / eff_Sel) ** 2 + (err_Glo / eff_Glo) ** 2)
    return eff_tot, err_tot


def eff_ztomumu_glo(nHLT_, nSel_, nGlo_, nSta_, nTrk_):
    #  facctorised Z efficiency from tag and probe at global muon level
    eff_Glo, err_Glo = muTagAndProbeEff(nHLT_, nSel_ + nGlo_, nTrk_ + nSta_)

    if eff_Glo <= 0:
        return 0., 0.

    eff_tot = eff_Glo ** 2
    err_tot = 2 * eff_Glo * err_Glo
    return eff_tot, err_tot


def zEff(nReco, nAcc, nBoth):
    """
    MC Z Effciciency
    :param nReco:
    :param nAcc:
    :param nBoth: intersect of nReco and nAcc for uncertainty estimation
    :return:
    """
    if nAcc < 100 or nReco < 50:
        return 0., 0.
    eff = nReco / nAcc
    x = nReco - nBoth
    y = nAcc - nBoth
    err = np.sqrt(x / (y + nBoth) ** 2
                  + ((x + nBoth) ** 2 * y) / (y + nBoth) ** 4
                  + ((x - y) ** 2 * nBoth) / (y - nBoth) ** 4
                  )
    return eff, err


def zMCEff1D(ZReco, ZAcc, bins, obs, region='inclusive', sel=''):
    """
    compute Z efficiency from MC and Z->mumu efficiency from tag and probe
    :param ZReco: dataframe with events satisfying the reco accepance
    :param ZAcc:  dataframe with events satisfying the gen level acceptance
    :param bins:
    :param obs:  oberservable
    :param region: string for naming reasons
    :param sel: selection on ZReco and ZAcc
    :return:
    """
    print(">>> make differential efficiencies for "+obs+" in "+region+" region")

    if sel != '':
        ZReco = ZReco.query(sel)
        ZAcc = ZAcc.query(sel)

    ZEff = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]
    ZEff_sel = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]
    ZEff_glo = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]

    MuEff_HLT_Sel = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Sel_Glo = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Glo_StaOrTrk = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Trk_Sta = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Sta_Trk = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]

    ZMuMuEff = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    ZMuMuEff_sel = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    ZMuMuEff_glo = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]

    for i in range(0, len(bins)-1):
        bin_low = bins[i]
        bin_high = bins[i+1]
        sel = '{0}>{1} & {0}<{2}'.format(obs, bin_low, bin_high)

        ZReco_i = ZReco.query(sel)

        qHLT = 'muon_Category==5 & antiMuon_Category==5'
        qSel = '(muon_Category==5 & antiMuon_Category==4) | (muon_Category==4 & antiMuon_Category==5)'
        qGlo = '(muon_Category==5 & antiMuon_Category==3) | (muon_Category==3 & antiMuon_Category==5)'
        qSta = '(muon_Category==5 & antiMuon_Category==2) | (muon_Category==2 & antiMuon_Category==5)'
        qTrk = '(muon_Category==5 & antiMuon_Category==1) | (muon_Category==1 & antiMuon_Category==5)'

        nHLT = len(ZReco_i.query(qHLT))
        nSel = len(ZReco_i.query(qSel))
        nGlo = len(ZReco_i.query(qGlo))
        nSta = len(ZReco_i.query(qSta))
        nTrk = len(ZReco_i.query(qTrk))

        MuEff_HLT_Sel[0][i], MuEff_HLT_Sel[1][i] = muTagAndProbeEff(nHLT, 0, nSel)
        MuEff_Sel_Glo[0][i], MuEff_Sel_Glo[1][i] = muTagAndProbeEff(nHLT, nSel, nGlo)
        MuEff_Glo_StaOrTrk[0][i], MuEff_Glo_StaOrTrk[1][i] = muTagAndProbeEff(nHLT, nSel + nGlo, nSta + nTrk)
        MuEff_Trk_Sta[0][i], MuEff_Trk_Sta[1][i] = muTagAndProbeEff(nHLT, nSel + nGlo + nTrk, nSta)
        MuEff_Sta_Trk[0][i], MuEff_Sta_Trk[1][i] = muTagAndProbeEff(nHLT, nSel + nGlo + nSta, nTrk)

        ZMuMuEff[0][i], ZMuMuEff[1][i] = eff_ztomumu_old(nHLT, nSel, nGlo, nSta, nTrk)
        ZMuMuEff_sel[0][i], ZMuMuEff_sel[1][i] = eff_ztomumu_sel(nHLT, nSel, nGlo, nSta, nTrk)
        ZMuMuEff_glo[0][i], ZMuMuEff_glo[1][i] = eff_ztomumu_glo(nHLT, nSel, nGlo, nSta, nTrk)

        #  full true Z efficiency
        ZAcc_i = ZAcc.query(sel)
        ZReco_i = ZReco.query(sel + '& (muon_Category==5 | antiMuon_Category==5) & muon_Category>=4 & antiMuon_Category>=4')

        nReco_i = len(ZReco_i)
        nAcc_i = len(ZAcc_i)
        nBoth_i = len(np.intersect1d(ZReco_i['LepPt'], ZAcc_i['LepPt']))

        ZEff[0][i], ZEff[1][i] = zEff(nReco_i, nAcc_i, nBoth_i)

        #  true Z efficiency at selection level (without requiring trigger)
        ZAcc_i = ZAcc.query(sel)
        ZReco_i = ZReco.query(sel + ' & muon_Category>=4 & antiMuon_Category>=4')

        nReco_i = len(ZReco_i)
        nAcc_i = len(ZAcc_i)
        nBoth_i = len(np.intersect1d(ZReco_i['LepPt'], ZAcc_i['LepPt']))

        ZEff_sel[0][i], ZEff_sel[1][i] = zEff(nReco_i, nAcc_i, nBoth_i)

        #  true Z efficiency at gobal muon level (requiring two global muons)
        ZAcc_i = ZAcc.query(sel)
        ZReco_i = ZReco.query(sel + '& muon_Category>=3 & antiMuon_Category>=3 ')

        nReco_i = len(ZReco_i)
        nAcc_i = len(ZAcc_i)
        nBoth_i = len(np.intersect1d(ZReco_i['LepPt'], ZAcc_i['LepPt']))

        ZEff_glo[0][i], ZEff_glo[1][i] = zEff(nReco_i, nAcc_i, nBoth_i)

    x = bins[:-1] + (bins[1:] - bins[:-1])/2

    def plot_zeff(eff_true, eff_tnp, name):
        plt.clf()
        fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        fig.subplots_adjust(hspace=0)

        ax[0].errorbar(x, eff_true[0], xerr=np.zeros(len(x)), yerr=eff_true[1], fmt='bo', label='true')
        ax[0].errorbar(x, eff_tnp[0], xerr=np.zeros(len(x)), yerr=eff_tnp[1], fmt='ro', label='tnp')
        ax[0].set(xlim=(bins[0], bins[-1]), ylim=(0.75, 1.1))
        ax[0].legend()
        ax[0].text(0.05, 0.9, r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$', transform=ax[0].transAxes)
        ax[0].text(0.05, 0.82, r'$p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4$', transform=ax[0].transAxes)
        ax[0].text(0.05, 0.74, region, transform=ax[0].transAxes)
        ax[0].set_ylabel(r'$\epsilon_\mathrm{Z}$')
        ax[0].set_xlabel(obs)
        ax[0].set_yticks([0.8, 0.85, 0.9, 0.95, 1.0, 1.05])
        ax[0].set_title("Z efficiency at {0} muon level".format(name))

        pulls = eff_tnp[0] - eff_true[0]
        pulls_sig = np.sqrt(eff_tnp[1]**2 + eff_true[1]**2)

        ax[1].errorbar(x, pulls, xerr=np.zeros(len(x)), yerr=pulls_sig, fmt='ko')#  , label='factorized - true')
        ax[1].plot(x, np.zeros(len(x)), color='gray', linestyle='dashed')
        ax[1].set_ylim(-0.02, 0.04)
        ax[1].set_ylabel(r'$\epsilon^\mathrm{tnp}_\mathrm{Z} - \epsilon^\mathrm{true}_\mathrm{Z}$')
        ax[1].set_xlabel(obs)
        ax[1].set_yticks([-0.01, 0., 0.01, 0.02, 0.03])
        plt.savefig(output+'/ZMuMu_{0}_{1}_{2}_level.png'.format(obs, region, name))
        plt.close()

    plot_zeff(ZEff, ZMuMuEff, "hlt")
    plot_zeff(ZEff_sel, ZMuMuEff_sel, "selection")
    plot_zeff(ZEff_glo, ZMuMuEff_glo, "global")



# acceptance selection
selection_Reco = 'z_recoMass > 66 ' \
                 '& z_recoMass < 116 ' \
                 '& muon_recoPt > 27 ' \
                 '& antiMuon_recoPt > 27 ' \
                 '& abs(muon_recoEta) < 2.4 ' \
                 '& abs(antiMuon_recoEta) < 2.4 '

selection_Gen = 'z_genMass > 66 ' \
                '& z_genMass < 116 ' \
                '& muon_genPt > 27 ' \
                '& antiMuon_genPt > 27 ' \
                '& abs(muon_genEta) < 2.4 ' \
                '& abs(antiMuon_genEta) < 2.4 '

# specify which branches to load
branches_Reco = ['nPV',
                 'muon_Category', 'antiMuon_Category',
                 'muon_recoEta', 'antiMuon_recoEta',
                 'muon_recoPhi', 'antiMuon_recoPhi',
                 'muon_recoPt', 'antiMuon_recoPt'
                 ]

branches_Gen = ['nPV',
                'muon_genEta', 'antiMuon_genEta',
                'muon_genPhi', 'antiMuon_genPhi',
                'muon_genPt', 'antiMuon_genPt'
                ]

print(">>> Load Events in reco acceptance")
dfReco = [tree_to_df(root2array(i, treeName[0], selection=selection_Reco, branches=branches_Reco), 5) for i in inputs]
dfReco = pd.concat(dfReco)
dfReco = dfReco.rename(columns={"muon_recoEta": "LepEta", "muon_recoPhi": "LepPhi", "muon_recoPt": "LepPt",
                                "antiMuon_recoEta": "AntiLepEta", "antiMuon_recoPhi": "AntiLepPhi", "antiMuon_recoPt": "AntiLepPt"
                                })

print(">>> Load Events in gen acceptance")
dfGen = [tree_to_df(root2array(i, treeName[0], selection=selection_Gen, branches=branches_Gen), 5) for i in inputs]
dfGen = pd.concat(dfGen)
dfGen = dfGen.rename(columns={"muon_genEta": "LepEta", "muon_genPhi": "LepPhi", "muon_genPt": "LepPt",
                              "antiMuon_genEta": "AntiLepEta", "antiMuon_genPhi": "AntiLepPhi", "antiMuon_genPt": "AntiLepPt"
                              })

print(">>> add new columns")
for df in (dfReco, dfGen):
    df['delPhiLL'] = abs(abs(df['LepPhi'] - df['AntiLepPhi']).apply(lambda x: x - 2*math.pi if x > math.pi else x))
    df['delEtaLL'] = abs(df['LepEta'] - df['AntiLepEta'])
    df['delRLL'] = np.sqrt((df['LepEta'] - df['AntiLepEta'])**2 + (df['LepPhi'] - df['AntiLepPhi'])**2)
    df['delPtLL'] = abs(df['LepPt'] - df['AntiLepPt'])
    df['relPtLL'] = abs(df['LepPt'] - df['AntiLepPt']) / abs(df['LepPt'] + df['AntiLepPt'])
    df['sumPtLL'] = df['LepPt'] + df['AntiLepPt']


# --- differential efficiencies
dfList = []
zMCEff1D(dfReco, dfGen, np.linspace(0., 5., 20), 'delRLL')
zMCEff1D(dfReco, dfGen, np.linspace(0., 5., 20), 'delRLL', 'BB', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 5., 20), 'delRLL', 'BE', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) > 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 5., 20), 'delRLL', 'EB', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 5., 20), 'delRLL', 'EE', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) > 0.9')

zMCEff1D(dfReco, dfGen, np.linspace(0., 2.5, 20), 'delEtaLL')
zMCEff1D(dfReco, dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'BB', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'BE', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) > 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'EB', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 2.5, 20), 'delEtaLL', 'EE', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) > 0.9')

zMCEff1D(dfReco, dfGen, np.linspace(0., math.pi, 20), 'delPhiLL')
zMCEff1D(dfReco, dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'BB', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'BE', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) > 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'EB', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., math.pi, 20), 'delPhiLL', 'EE', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) > 0.9')

zMCEff1D(dfReco, dfGen, np.linspace(0., 100, 20), 'delPtLL')
zMCEff1D(dfReco, dfGen, np.linspace(0., 100, 20), 'delPtLL', 'BB', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 100, 20), 'delPtLL', 'BE', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) > 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 100, 20), 'delPtLL', 'EB', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 100, 20), 'delPtLL', 'EE', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) > 0.9')

zMCEff1D(dfReco, dfGen, np.linspace(0., 1., 20), 'relPtLL')
zMCEff1D(dfReco, dfGen, np.linspace(0., 1., 20), 'relPtLL', 'BB', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 1., 20), 'relPtLL', 'BE', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) > 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 1., 20), 'relPtLL', 'EB', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0., 1., 20), 'relPtLL', 'EE', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) > 0.9')

zMCEff1D(dfReco, dfGen, np.linspace(50., 250, 20), 'sumPtLL')
zMCEff1D(dfReco, dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'BB', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'BE', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) > 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'EB', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(50., 250, 20), 'sumPtLL', 'EE', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) > 0.9')

zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'inclusive')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'BB', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'BE', sel='abs(LepEta) < 0.9 & abs(AntiLepEta) > 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'EB', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) < 0.9')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'EE', sel='abs(LepEta) > 0.9 & abs(AntiLepEta) > 0.9')

zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delRLL0to25', sel='delRLL < 2.5')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delRLL25toInf', sel='delRLL > 2.5')

zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delEtaLL0to1', sel='delEtaLL < 1')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delEtaLL1toInf', sel='delEtaLL > 1')

zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPhiLL0to25', sel='delPhiLL < 2.5')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPhiLL25toPi', sel='delPhiLL > 2.5')

zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPtLL0to10', sel='delPtLL < 15')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'delPtLL10toInf', sel='delPtLL > 15')

zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'relPtLL0to01', sel='relPtLL < 0.1')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'relPtLL01toInf', sel='relPtLL > 0.1')

zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'sumPtLL0to100', sel='sumPtLL < 100')
zMCEff1D(dfReco, dfGen, np.linspace(0.5, 59.5, 60), 'nPV', 'sumPtLL100toInf', sel='sumPtLL > 100')
