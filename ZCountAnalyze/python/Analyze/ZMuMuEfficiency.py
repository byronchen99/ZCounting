from __future__ import division, print_function
# # #
#
#  Measures the Z->mumu efficiency by factorization from single muon efficiencies
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
selection='ZStableMass > 66 ' \
          '& ZStableMass < 116 ' \
          '& (ZDecayMode == 13 | ZDecayMode == 151313) ' \
          '& ZLeptonPt > 27 ' \
          '& ZAntiLeptonPt > 27 ' \
          '& abs(ZLeptonEta) < 2.4 ' \
          '& abs(ZAntiLeptonEta) < 2.4 '
#specify which branches to load
branches=['MuonProbeCategory','nPV','eventWeight']

from Utils import tree_to_df
ZAcc = tree_to_df(root2array(input, treeName[0], selection=selection, branches=branches),5)
ZAcc['eventWeight'] = ZAcc['eventWeight']/abs(ZAcc['eventWeight'])

### alysis part ###

#dataframe with all Zs that are reconstructed
#ZReco = ZAcc.query('MuonProbeCategory_0==1 | MuonProbeCategory_0==2')

nHLT = np.sum(ZAcc.query('MuonProbeCategory_0==1')['eventWeight'])
nSel = np.sum(ZAcc.query('MuonProbeCategory_0==2')['eventWeight'])
nGlo = np.sum(ZAcc.query('MuonProbeCategory_0==3')['eventWeight'])
nSta = np.sum(ZAcc.query('MuonProbeCategory_0==4')['eventWeight'])
nTrk = np.sum(ZAcc.query('MuonProbeCategory_0==5')['eventWeight'])

pdb.set_trace()

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

#Efficiency in dependence of PU
minPU=0
maxPU=60

ZEff_PU_old = np.zeros(maxPU-minPU)
ZEff_PU_old_err = np.zeros(maxPU-minPU)
ZEff_PU_new = np.zeros(maxPU-minPU)
ZEff_PU_new_err = np.zeros(maxPU-minPU)

for iPV in range(minPU,maxPU):
    nHLT = np.sum(ZAcc.query('MuonProbeCategory_0==1 & nPV=={0}'.format(iPV))['eventWeight'])
    nSel = np.sum(ZAcc.query('MuonProbeCategory_0==2 & nPV=={0}'.format(iPV))['eventWeight'])
    nGlo = np.sum(ZAcc.query('MuonProbeCategory_0==3 & nPV=={0}'.format(iPV))['eventWeight'])
    nSta = np.sum(ZAcc.query('MuonProbeCategory_0==4 & nPV=={0}'.format(iPV))['eventWeight'])
    nTrk = np.sum(ZAcc.query('MuonProbeCategory_0==5 & nPV=={0}'.format(iPV))['eventWeight'])

    ZEff_PU_old[iPV], ZEff_PU_old_err[iPV] = EffMuMu_old(nHLT,nSel,nGlo,nSta,nTrk)
    ZEff_PU_new[iPV], ZEff_PU_new_err[iPV] = EffMuMu_new(nHLT,nSel,nGlo,nSta,nTrk)


dfOut = DataFrame()
dfOut['nPV'] = range(minPU,maxPU)
dfOut['ZMuMuEff_old'] = ZEff_PU_old
dfOut['ZMuMuErr_old'] = ZEff_PU_old_err
dfOut['ZMuMuEff_new'] = ZEff_PU_new
dfOut['ZMuMuErr_new'] = ZEff_PU_new_err
dfOut.to_hdf('ZMuMuEfficiencies.h5', 'ZMuMuEfficiencies', mode='w')


### Plot ###

plt.errorbar(range(minPU,maxPU), ZEff_PU_old, xerr=np.zeros(maxPU-minPU), yerr=ZEff_PU_old_err, fmt='ro' )
plt.ylim((0.7,1.1))
plt.ylabel('old Z efficiency')
plt.xlabel('nPV')
plt.savefig('ZMuMuEff_PU_old.pdf')
plt.clf()

plt.errorbar(range(minPU,maxPU), ZEff_PU_new, xerr=np.zeros(maxPU-minPU), yerr=ZEff_PU_new_err, fmt='go' )
plt.ylim((0.7,1.1))
plt.ylabel('new Z efficiency')
plt.xlabel('nPV')
plt.savefig('ZMuMuEff_PU_new.pdf')
plt.clf()

plt.errorbar(range(minPU,maxPU), ZEff_PU_old, xerr=np.zeros(maxPU-minPU), yerr=ZEff_PU_old_err, fmt='ro', label='old' )
plt.errorbar(range(minPU,maxPU), ZEff_PU_new, xerr=np.zeros(maxPU-minPU), yerr=ZEff_PU_new_err, fmt='go', label='new' )
plt.ylim((0.7,1.1))
plt.legend()
plt.ylabel('Z efficiency')
plt.xlabel('nPV')
plt.savefig('ZMuMuEff_PU_both.pdf')