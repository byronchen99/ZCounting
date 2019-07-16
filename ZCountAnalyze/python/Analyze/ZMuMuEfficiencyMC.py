from __future__ import division, print_function
# # #
#
#  Measures the Z->mumu efficiency by factorization from single muon efficiencies
#    using gen level matched reco particles
#
# # #
from root_numpy import root2array, list_trees
import pandas as pd
import pdb
import argparse
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='./ZMuMUEfficiency')
parser.add_argument(
    '-i','--input', nargs='+',
    help='specify input root file'
)
args = parser.parse_args()

inputs = args.input

if isinstance(inputs, (list,)):
    treeName = list_trees(inputs[0])
else:
    treeName = list_trees(inputs)
    inputs = [inputs,]

if(len(treeName) > 1):
    print("more then one tree in file ... specify, which tree to use")
    exit()

#acceptance selection
selection='ZMassReco > 66 ' \
          '& ZMassReco < 116 ' \
          '& (ZDecayMode == 13 | ZDecayMode == 151313) ' \
          '& ZLeptonRecoPt > 27 ' \
          '& ZAntiLeptonRecoPt > 27 ' \
          '& abs(ZLeptonRecoEta) < 2.4 ' \
          '& abs(ZAntiLeptonRecoEta) < 2.4 '
#specify which branches to load
branches=['ZLeptonRecoCat','ZAntiLeptonRecoCat','nPV','eventWeight', 'ZLeptonRecoEta', 'ZAntiLeptonRecoEta']

from Utils import tree_to_df
ZAcc = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches),5) for i in inputs]
ZAcc = pd.concat(ZAcc)
ZAcc['eventWeight'] = ZAcc['eventWeight']/abs(ZAcc['eventWeight'])

### alysis part ###

#dataframe with all Zs that are reconstructed
#ZReco = ZAcc.query('MuonProbeCategory_0==1 | MuonProbeCategory_0==2')

qHLT = 'ZLeptonRecoCat==1 & ZAntiLeptonRecoCat==1'
qSel = '(ZLeptonRecoCat==1 & ZAntiLeptonRecoCat == 2) | (ZLeptonRecoCat==2 & ZAntiLeptonRecoCat == 1)'
qGlo = '(ZLeptonRecoCat==1 & ZAntiLeptonRecoCat == 3) | (ZLeptonRecoCat==3 & ZAntiLeptonRecoCat == 1)'
qSta = '(ZLeptonRecoCat==1 & ZAntiLeptonRecoCat == 4) | (ZLeptonRecoCat==4 & ZAntiLeptonRecoCat == 1)'
qTrk = '(ZLeptonRecoCat==1 & ZAntiLeptonRecoCat == 5) | (ZLeptonRecoCat==5 & ZAntiLeptonRecoCat == 1)'

nHLT = np.sum(ZAcc.query(qHLT)['eventWeight'])
nSel = np.sum(ZAcc.query(qSel)['eventWeight'])
nGlo = np.sum(ZAcc.query(qGlo)['eventWeight'])
nSta = np.sum(ZAcc.query(qSta)['eventWeight'])
nTrk = np.sum(ZAcc.query(qTrk)['eventWeight'])

def Eff(hlt,pp, fp):
    if 2*hlt+pp <= 0:
        return 0., 0.
    elif fp <= 0:
        return 0., 0.

    eff = (2*hlt+pp)/(2*hlt+pp+fp)
    eff_err = 1. / (2*hlt+fp+pp)**2 * np.sqrt( pp*fp**2 + hlt*(2*fp)**2 + fp*(pp+2*hlt)**2 )
    return eff, eff_err

def EffMuMu_old(nHLT_,nSel_,nGlo_,nSta_,nTrk_):
    eff_HLT, err_HLT = Eff(nHLT_, 0, nSel_)
    eff_Sel, err_Sel = Eff(nHLT_, nSel_, nGlo_)
    eff_Glo, err_Glo = Eff(nHLT_, nSel_ + nGlo_, nTrk_ + nSta_)

    if eff_HLT <= 0 or eff_Sel <=0 or eff_Glo <= 0:
        return 0., 0.

    eff_tot = (1 - (1 - eff_HLT)**2) * eff_Sel**2 * eff_Glo**2
    err_tot = 2 * eff_tot * np.sqrt(((1-eff_HLT)/(1-(1-eff_HLT)**2)*err_HLT)**2 + (err_Sel/eff_Sel)**2 + (err_Glo/eff_Glo)**2)
    return eff_tot, err_tot


def EffMuMu_new(nHLT_, nSel_, nGlo_, nSta_, nTrk_):
    eff_HLT, err_HLT = Eff(nHLT_, 0, nSel_)
    eff_Sel, err_Sel = Eff(nHLT_, nSel_, nGlo_)
    eff_Trk, err_Trk = Eff(nHLT_, nSel_ + nGlo_ + nTrk_, nSta_)
    eff_Sta, err_Sta = Eff(nHLT_, nSel_ + nGlo_ + nSta_, nTrk_)

    if eff_HLT <= 0 or eff_Sel <=0 or eff_Trk <= 0 or eff_Sta <= 0:
        return 0., 0.

    eff_tot = (1 - (1 - eff_HLT)**2) * eff_Sel**2 * eff_Trk**2 * eff_Sta**2
    err_tot = 2 * eff_tot * np.sqrt(((1-eff_HLT)/(1-(1-eff_HLT)**2)*err_HLT)**2 + (err_Sel/eff_Sel)**2 + (err_Trk/eff_Trk)**2 + (err_Sta/eff_Sta)**2)
    return eff_tot, err_tot


#old Z->MuMu Efficiency
print("inclusive Z->mumu old Efficiency is: eff = "+str(EffMuMu_old(nHLT,nSel,nGlo,nSta,nTrk)))

#new Z->MuMu Efficiency with tracking
print("inclusive Z->mumu new Efficiency is: eff = "+str(EffMuMu_new(nHLT,nSel,nGlo,nSta,nTrk)))

### Plot ###
def plotEff(name,xlabel,ylabel,bins,eff,err, cuts):
    plt.clf()
    fig, ax = plt.subplots()

    ax.errorbar(bins[:-1] + (bins[1:] - bins[:-1])/2 , eff, xerr=np.zeros(len(bins)-1), yerr=err, fmt='ro')
    ax.set(xlim=(bins[0], bins[-1]), ylim=(0.7, 1.1))
    ax.text(0.05, 0.9, cuts, transform=ax.transAxes)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    plt.savefig(name+'.png')
    plt.close()


#Differential efficiencies

dfOut = pd.DataFrame()


def differentialEff(bins, obs, cuts, sel=''):
    if sel != '':
        sel = '&'+sel
    
    MuEff_HLT_Sel = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]
    MuEff_Sel_Glo = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]
    MuEff_Glo_StaOrTrk = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]
    MuEff_Trk_Sta = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]
    MuEff_Sta_Trk = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]

    ZEff_old = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]
    ZEff_new = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]

    for i in range(0,len(bins)-1):
        bin_low = bins[i]
        bin_high = bins[i+1]
        nHLT = np.sum(ZAcc.query('('+qHLT+') & {0}>{1} & {0}<{2} {3}'.format(obs, bin_low, bin_high, sel))['eventWeight'])
        nSel = np.sum(ZAcc.query('('+qSel+') & {0}>{1} & {0}<{2} {3}'.format(obs, bin_low, bin_high, sel))['eventWeight'])
        nGlo = np.sum(ZAcc.query('('+qGlo+') & {0}>{1} & {0}<{2} {3}'.format(obs, bin_low, bin_high, sel))['eventWeight'])
        nSta = np.sum(ZAcc.query('('+qSta+') & {0}>{1} & {0}<{2} {3}'.format(obs, bin_low, bin_high, sel))['eventWeight'])
        nTrk = np.sum(ZAcc.query('('+qTrk+') & {0}>{1} & {0}<{2} {3}'.format(obs, bin_low, bin_high, sel))['eventWeight'])
    
        MuEff_HLT_Sel[0][i], MuEff_HLT_Sel[1][i] = Eff(nHLT, 0, nSel)
        MuEff_Sel_Glo[0][i], MuEff_Sel_Glo[1][i] = Eff(nHLT, nSel, nGlo)
        MuEff_Glo_StaOrTrk[0][i], MuEff_Glo_StaOrTrk[1][i] = Eff(nHLT, nSel + nGlo, nSta + nTrk)
        MuEff_Trk_Sta[0][i], MuEff_Trk_Sta[1][i] = Eff(nHLT, nSel + nGlo + nTrk, nSta)
        MuEff_Sta_Trk[0][i], MuEff_Sta_Trk[1][i] = Eff(nHLT, nSel + nGlo + nSta, nTrk)
    
        ZEff_old[0][i], ZEff_old[1][i] = EffMuMu_old(nHLT,nSel,nGlo,nSta,nTrk)
        ZEff_new[0][i], ZEff_new[1][i] = EffMuMu_new(nHLT,nSel,nGlo,nSta,nTrk)

    dfOut[obs+'_'+cuts+'_lowEdge'] = bins[:-1]
    dfOut[obs+'_'+cuts+'_upEdge'] = bins[1:]

    dfOut['ZMuMuEff_{0}_{1}_old'.format(obs, cuts)] = ZEff_old[0]
    dfOut['ZMuMuErr_{0}_{1}_old'.format(obs, cuts)] = ZEff_old[1]
    dfOut['ZMuMuEff_{0}_{1}_new'.format(obs, cuts)] = ZEff_new[0]
    dfOut['ZMuMuErr_{0}_{1}_new'.format(obs, cuts)] = ZEff_new[1]

    plotEff('MuEff_{0}_{1}_HLT_Sel'.format(obs, cuts), obs, 'Mu HLT Efficiency',bins, MuEff_HLT_Sel[0], MuEff_HLT_Sel[1], cuts)
    plotEff('MuEff_{0}_{1}_Sel_Glo'.format(obs, cuts), obs, 'Mu Sel Efficiency',bins, MuEff_Sel_Glo[0], MuEff_Sel_Glo[1], cuts)
    plotEff('MuEff_{0}_{1}_Glo_StaOrTrk'.format(obs, cuts), obs, 'Mu Glo Efficiency',bins, MuEff_Glo_StaOrTrk[0], MuEff_Glo_StaOrTrk[1], cuts)
    plotEff('MuEff_{0}_{1}_Trk_Sta'.format(obs, cuts), obs, 'Mu Trk Efficiency',bins, MuEff_Trk_Sta[0], MuEff_Trk_Sta[1], cuts)
    plotEff('MuEff_{0}_{1}_Sta_Trk'.format(obs, cuts), obs, 'Mu Sta Efficiency',bins, MuEff_Sta_Trk[0], MuEff_Sta_Trk[1], cuts)

    plotEff('ZEff_{0}_{1}_old'.format(obs, cuts), obs, 'old Z efficiency',bins, ZEff_old[0], ZEff_old[1], cuts)
    plotEff('ZEff_{0}_{1}_new'.format(obs, cuts), obs, 'new Z efficiency',bins, ZEff_new[0], ZEff_new[1], cuts)

    plt.clf()
    fig, ax = plt.subplots()
    #fig.subplots_adjust(bottom=0.15, left=0.2)

    ax.errorbar(bins[:-1] + (bins[1:] - bins[:-1])/2, ZEff_old[0], xerr=np.zeros(len(bins)-1), yerr=ZEff_old[1], fmt='ro', label='old')
    ax.errorbar(bins[:-1] + (bins[1:] - bins[:-1])/2, ZEff_new[0], xerr=np.zeros(len(bins)-1), yerr=ZEff_new[1], fmt='go', label='new')
    ax.set(xlim=(bins[0], bins[-1]), ylim=(0.7, 1.1))
    ax.legend()
    ax.text(0.05, 0.9, cuts, transform=ax.transAxes)
    ax.set_ylabel('Z efficiency')
    ax.set_xlabel(obs)
    plt.savefig('ZMuMu_{0}_{1}_both.png'.format(obs,cuts))
    plt.close()


differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'inclusive')
differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'BB', sel='abs(ZLeptonRecoEta) < 0.9 & abs(ZAntiLeptonRecoEta) < 0.9')
differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'EB', sel='abs(ZLeptonRecoEta) > 0.9 & abs(ZAntiLeptonRecoEta) < 0.9')
differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'BE', sel='abs(ZLeptonRecoEta) < 0.9 & abs(ZAntiLeptonRecoEta) > 0.9')
differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'EE', sel='abs(ZLeptonRecoEta) > 0.9 & abs(ZAntiLeptonRecoEta) > 0.9')

dfOut.to_hdf('ZMuMuEfficiencies.h5', 'ZMuMuEfficiencies', mode='w')

