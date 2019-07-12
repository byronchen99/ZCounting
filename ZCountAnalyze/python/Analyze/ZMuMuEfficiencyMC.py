from __future__ import division, print_function
# # #
#
#  Measures the Z->mumu efficiency by factorization from single muon efficiencies
#    using gen level matched reco particles
#
# # #
from root_numpy import root2array, list_trees
from pandas import DataFrame
import pdb
import argparse
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(prog='./ZMuMUEfficiency')
parser.add_argument(
    '-i','--input', required=True,
    help='specify input root file'
)
args = parser.parse_args()

input = args.input

treeName = list_trees(input)
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
branches=['ZLeptonRecoCat','ZAntiLeptonRecoCat','nPV','eventWeight']

from Utils import tree_to_df
ZAcc = tree_to_df(root2array(input, treeName[0], selection=selection, branches=branches),5)
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

#Differential efficiencies
minPU=0
maxPU=60

MuEff_PU_HLT_Sel = [np.zeros(maxPU-minPU), np.zeros(maxPU-minPU)]
MuEff_PU_Sel_Glo = [np.zeros(maxPU-minPU), np.zeros(maxPU-minPU)]
MuEff_PU_Glo_StaOrTrk = [np.zeros(maxPU-minPU), np.zeros(maxPU-minPU)]
MuEff_PU_Trk_Sta = [np.zeros(maxPU-minPU), np.zeros(maxPU-minPU)]
MuEff_PU_Sta_Trk = [np.zeros(maxPU-minPU), np.zeros(maxPU-minPU)]

ZEff_PU_old = [np.zeros(maxPU-minPU), np.zeros(maxPU-minPU)]
ZEff_PU_new = [np.zeros(maxPU-minPU), np.zeros(maxPU-minPU)]


for iPV in range(minPU,maxPU):
    nHLT = np.sum(ZAcc.query(qHLT+' & nPV=={0}'.format(iPV))['eventWeight'])
    nSel = np.sum(ZAcc.query(qSel+' & nPV=={0}'.format(iPV))['eventWeight'])
    nGlo = np.sum(ZAcc.query(qGlo+' & nPV=={0}'.format(iPV))['eventWeight'])
    nSta = np.sum(ZAcc.query(qSta+' & nPV=={0}'.format(iPV))['eventWeight'])
    nTrk = np.sum(ZAcc.query(qTrk+' & nPV=={0}'.format(iPV))['eventWeight'])

    MuEff_PU_HLT_Sel[0][iPV], MuEff_PU_HLT_Sel[1][iPV] = Eff(nHLT, 0, nSel)
    MuEff_PU_Sel_Glo[0][iPV], MuEff_PU_Sel_Glo[1][iPV] = Eff(nHLT, nSel, nGlo)
    MuEff_PU_Glo_StaOrTrk[0][iPV], MuEff_PU_Glo_StaOrTrk[1][iPV] = Eff(nHLT, nSel + nGlo, nSta + nTrk)
    MuEff_PU_Trk_Sta[0][iPV], MuEff_PU_Trk_Sta[1][iPV] = Eff(nHLT, nSel + nGlo + nTrk, nSta)
    MuEff_PU_Sta_Trk[0][iPV], MuEff_PU_Sta_Trk[1][iPV] = Eff(nHLT, nSel + nGlo + nSta, nTrk)

    ZEff_PU_old[0][iPV], ZEff_PU_old[1][iPV] = EffMuMu_old(nHLT,nSel,nGlo,nSta,nTrk)
    ZEff_PU_new[0][iPV], ZEff_PU_new[1][iPV] = EffMuMu_new(nHLT,nSel,nGlo,nSta,nTrk)


dfOut = DataFrame()
dfOut['nPV'] = range(minPU,maxPU)
dfOut['ZMuMuEff_old'] = ZEff_PU_old[0]
dfOut['ZMuMuErr_old'] = ZEff_PU_old[1]
dfOut['ZMuMuEff_new'] = ZEff_PU_new[0]
dfOut['ZMuMuErr_new'] = ZEff_PU_new[1]
dfOut.to_hdf('ZMuMuEfficiencies.h5', 'ZMuMuEfficiencies', mode='w')


### Plot ###
def plotEff(name,xlabel,ylabel,xMin,xMax,eff,err=None):
    plt.clf()
    if not isinstance(err, (list,)):
        err = np.zeros(xMax - xMin)

    plt.errorbar(range(xMin, xMax), eff, xerr=np.zeros(xMax - xMin), yerr=err, fmt='ro')
    plt.ylim((0.7, 1.1))
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.savefig(name+'.pdf')

plotEff('MuEff_PU_HLT_Sel','nPV', 'Mu HLT Efficiency',minPU, maxPU, MuEff_PU_HLT_Sel[0], MuEff_PU_HLT_Sel[1])
plotEff('MuEff_PU_Sel_Glo','nPV', 'Mu Sel Efficiency',minPU, maxPU, MuEff_PU_Sel_Glo[0], MuEff_PU_Sel_Glo[1])
plotEff('MuEff_PU_Glo_StaOrTrk','nPV', 'Mu Glo Efficiency',minPU, maxPU, MuEff_PU_Glo_StaOrTrk[0], MuEff_PU_Glo_StaOrTrk[1])
plotEff('MuEff_PU_Trk_Sta','nPV', 'Mu Trk Efficiency',minPU, maxPU, MuEff_PU_Trk_Sta[0], MuEff_PU_Trk_Sta[1])
plotEff('MuEff_PU_Sta_Trk','nPV', 'Mu Sta Efficiency',minPU, maxPU, MuEff_PU_Sta_Trk[0], MuEff_PU_Sta_Trk[1])

plotEff('ZEff_PU_old','nPV', 'old Z efficiency',minPU, maxPU, ZEff_PU_old[0], ZEff_PU_old[1])
plotEff('ZEff_PU_new','nPV', 'new Z efficiency',minPU, maxPU, ZEff_PU_new[0], ZEff_PU_new[1])

plt.errorbar(range(minPU,maxPU), ZEff_PU_old[0], xerr=np.zeros(maxPU-minPU), yerr=ZEff_PU_old[1], fmt='ro', label='old' )
plt.errorbar(range(minPU,maxPU), ZEff_PU_new[0], xerr=np.zeros(maxPU-minPU), yerr=ZEff_PU_new[1], fmt='go', label='new' )
plt.ylim((0.7,1.1))
plt.legend()
plt.ylabel('Z efficiency')
plt.xlabel('nPV')
plt.savefig('ZMuMuEff_PU_both.pdf')