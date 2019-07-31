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


def ztoMuMuEff_old(nHLT_,nSel_,nGlo_,nSta_,nTrk_):
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


def ztoMuMuEff_new(nHLT_, nSel_, nGlo_, nSta_, nTrk_):
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

    MuEff_HLT_Sel = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Sel_Glo = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Glo_StaOrTrk = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Trk_Sta = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]
    MuEff_Sta_Trk = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]

    ZMuMuEff = [np.zeros(len(bins) - 1), np.zeros(len(bins) - 1)]

    for i in range(0, len(bins)-1):
        bin_low = bins[i]
        bin_high = bins[i+1]
        sel = '{0}>{1} & {0}<{2}'.format(obs, bin_low, bin_high)

        ZReco_i = ZReco.query(sel)

        qHLT = 'ZLeptonRecoCat_1==5 & ZAntiLeptonRecoCat_1==5'
        qSel = '(ZLeptonRecoCat_1==5 & ZAntiLeptonRecoCat_1 == 4) | (ZLeptonRecoCat_1==4 & ZAntiLeptonRecoCat_1 == 5)'
        qGlo = '(ZLeptonRecoCat_1==5 & ZAntiLeptonRecoCat_1 == 3) | (ZLeptonRecoCat_1==3 & ZAntiLeptonRecoCat_1 == 5)'
        qSta = '(ZLeptonRecoCat_1==5 & ZAntiLeptonRecoCat_1 == 2) | (ZLeptonRecoCat_1==2 & ZAntiLeptonRecoCat_1 == 5)'
        qTrk = '(ZLeptonRecoCat_1==5 & ZAntiLeptonRecoCat_1 == 1) | (ZLeptonRecoCat_1==1 & ZAntiLeptonRecoCat_1 == 5)'

        nHLT = np.sum(ZReco_i.query(qHLT)['eventWeight'])
        nSel = np.sum(ZReco_i.query(qSel)['eventWeight'])
        nGlo = np.sum(ZReco_i.query(qGlo)['eventWeight'])
        nSta = np.sum(ZReco_i.query(qSta)['eventWeight'])
        nTrk = np.sum(ZReco_i.query(qTrk)['eventWeight'])

        MuEff_HLT_Sel[0][i], MuEff_HLT_Sel[1][i] = muTagAndProbeEff(nHLT, 0, nSel)
        MuEff_Sel_Glo[0][i], MuEff_Sel_Glo[1][i] = muTagAndProbeEff(nHLT, nSel, nGlo)
        MuEff_Glo_StaOrTrk[0][i], MuEff_Glo_StaOrTrk[1][i] = muTagAndProbeEff(nHLT, nSel + nGlo, nSta + nTrk)
        MuEff_Trk_Sta[0][i], MuEff_Trk_Sta[1][i] = muTagAndProbeEff(nHLT, nSel + nGlo + nTrk, nSta)
        MuEff_Sta_Trk[0][i], MuEff_Sta_Trk[1][i] = muTagAndProbeEff(nHLT, nSel + nGlo + nSta, nTrk)

        ZMuMuEff[0][i], ZMuMuEff[1][i] = ztoMuMuEff_old(nHLT, nSel, nGlo, nSta, nTrk)

        ZAcc_i = ZAcc.query(sel)
        ZReco_i = ZReco_i.query('(ZAntiLeptonRecoCat_0>=4 & ZLeptonRecoCat_0==5) '
                                '| (ZAntiLeptonRecoCat_0==5 & ZLeptonRecoCat_0>=4)')

        nReco_i = np.sum(ZReco_i['eventWeight'])
        nAcc_i = np.sum(ZAcc_i['eventWeight'])
        nBoth_i = len(np.intersect1d(ZReco_i['LepPt'], ZAcc_i['LepPt']))

        ZEff[0][i], ZEff[1][i] = zEff(nReco_i, nAcc_i, nBoth_i)

    x = bins[:-1] + (bins[1:] - bins[:-1])/2

    plt.clf()
    fig, ax = plt.subplots()
    ax.errorbar(x, ZEff[0], xerr=np.zeros(len(x)), yerr=ZEff[1], fmt='bo', label='true')
    ax.errorbar(x, ZMuMuEff[0], xerr=np.zeros(len(x)), yerr=ZMuMuEff[1], fmt='ro', label='factorized')
    ax.set(xlim=(bins[0], bins[-1]), ylim=(0.8, 1.1))
    ax.legend()
    ax.text(0.05, 0.9, r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$', transform=ax.transAxes)
    ax.text(0.05, 0.85, r'$p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4$', transform=ax.transAxes)
    ax.text(0.05, 0.8, region, transform=ax.transAxes)
    ax.set_ylabel('Z efficiency')
    ax.set_xlabel(obs)
    plt.savefig(output+'/ZMuMu_{0}_{1}_all.png'.format(obs, region))
    plt.close()

    pulls = ZMuMuEff[0] - ZEff[0]
    pulls_sig = np.sqrt(ZMuMuEff[1]**2 + ZEff[1]**2)

    plt.clf()
    fig, ax = plt.subplots()
    ax.errorbar(x, pulls, xerr=np.zeros(len(x)), yerr=pulls_sig, fmt='ro', label='factorized - true')
    ax.set(xlim=(bins[0], bins[-1]), ylim=(-0.025, 0.05))
    ax.text(0.05, 0.9, r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$', transform=ax.transAxes)
    ax.text(0.05, 0.85, r'$p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4$', transform=ax.transAxes)
    ax.text(0.05, 0.8, region, transform=ax.transAxes)
    ax.legend(loc='upper right')
    ax.set_ylabel('Z efficiency pulls')
    ax.set_xlabel(obs)
    plt.savefig(output+'/ZMuMu_{0}_{1}_pulls.png'.format(obs, region))
    plt.close()


#    dfOut = pd.DataFrame()
#    dfOut[obs + '_' + region + '_lowEdge'] = bins[:-1]
#    dfOut[obs + '_' + region + '_upEdge'] = bins[1:]
#    dfOut['ZEff_{0}_{1}'.format(obs, region)] = ZEff[0]
#    dfOut['ZErr_{0}_{1}'.format(obs, region)] = ZEff[1]
#    dfOut['ZMuMuEff_{0}_{1}'.format(obs, region)] = ZMuMuEff[0]
#    dfOut['ZMuMuErr_{0}_{1}'.format(obs, region)] = ZMuMuEff[1]

#    return dfOut


# acceptance selection
selection_Reco = 'ZMassReco > 66 ' \
                 '& ZMassReco < 116 ' \
                 '& ZDecayMode == 13 ' \
                 '& ZLeptonRecoPt > 27 ' \
                 '& ZAntiLeptonRecoPt > 27 ' \
                 '& abs(ZLeptonRecoEta) < 2.4 ' \
                 '& abs(ZAntiLeptonRecoEta) < 2.4 '

selection_Gen = 'ZStableMass > 66 ' \
                '& ZStableMass < 116 ' \
                '& ZDecayMode == 13 ' \
                '& ZLeptonPt > 27 ' \
                '& ZAntiLeptonPt > 27 ' \
                '& abs(ZLeptonEta) < 2.4 ' \
                '& abs(ZAntiLeptonEta) < 2.4 '

# specify which branches to load
branches_Reco = ['ZLeptonRecoCat', 'ZAntiLeptonRecoCat', 'nPV', 'eventWeight',
                 'ZLeptonRecoEta', 'ZAntiLeptonRecoEta',
                 'ZLeptonRecoPhi', 'ZAntiLeptonRecoPhi',
                 'ZLeptonRecoPt', 'ZAntiLeptonRecoPt'
                 ]

branches_Gen = ['nPV', 'eventWeight',
                'ZLeptonEta', 'ZAntiLeptonEta',
                'ZLeptonPhi', 'ZAntiLeptonPhi',
                'ZLeptonPt', 'ZAntiLeptonPt'
                ]

print(">>> Load Events in reco acceptance")
dfReco = [tree_to_df(root2array(i, treeName[0], selection=selection_Reco, branches=branches_Reco), 5) for i in inputs]
dfReco = pd.concat(dfReco)
dfReco['eventWeight'] = dfReco['eventWeight']/dfReco['eventWeight']
dfReco = dfReco.rename(columns={"ZLeptonRecoEta": "LepEta", "ZLeptonRecoPhi": "LepPhi", "ZLeptonRecoPt": "LepPt",
                                "ZAntiLeptonRecoEta": "AntiLepEta", "ZAntiLeptonRecoPhi": "AntiLepPhi", "ZAntiLeptonRecoPt": "AntiLepPt"
                                })

print(">>> Load Events in gen acceptance")
dfGen = [tree_to_df(root2array(i, treeName[0], selection=selection_Gen, branches=branches_Gen), 5) for i in inputs]
dfGen = pd.concat(dfGen)
dfGen['eventWeight'] = dfGen['eventWeight']/dfGen['eventWeight']
dfGen = dfGen.rename(columns={"ZLeptonEta": "LepEta", "ZLeptonPhi": "LepPhi", "ZLeptonPt": "LepPt",
                              "ZAntiLeptonEta": "AntiLepEta", "ZAntiLeptonPhi": "AntiLepPhi", "ZAntiLeptonPt": "AntiLepPt"
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



