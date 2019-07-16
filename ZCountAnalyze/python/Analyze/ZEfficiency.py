from __future__ import division, print_function
# # #
#
#  Measures the true Z->mumu reconstruction efficiency from MC
#
# # #
from root_numpy import root2array, list_trees
import pdb
import argparse
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(prog='./ZEfficiency')
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
selection='ZStableMass > 66 ' \
          '& ZStableMass < 116 ' \
          '& (ZDecayMode == 13 | ZDecayMode == 151313) ' \
          '& ZLeptonPt > 27 ' \
          '& ZAntiLeptonPt > 27 ' \
          '& abs(ZLeptonEta) < 2.4 ' \
          '& abs(ZAntiLeptonEta) < 2.4 '
#specify which branches to load
branches=['nPV','eventWeight','ZLeptonRecoEta','ZAntiLeptonRecoEta','ZLeptonEta','ZAntiLeptonEta']

from Utils import tree_to_df

#dataframe with all Zs that are accepted in Gen
ZAcc = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches),5) for i in inputs]
ZAcc = pd.concat(ZAcc)
ZAcc['eventWeight'] = ZAcc['eventWeight']/abs(ZAcc['eventWeight'])
nAcc = np.sum(ZAcc['eventWeight'])

print("number of accepted gen Events is = "+str(nAcc))


#dataframe with all Zs that are reconstructed
selection='(' \
          '    (ZAntiLeptonRecoCat==1 & ZLeptonRecoCat==2) ' \
            '| (ZAntiLeptonRecoCat==2 & ZLeptonRecoCat==1) ' \
            '| (ZAntiLeptonRecoCat==1 & ZLeptonRecoCat==1)' \
          ') ' \
          '& ZMassReco > 66 ' \
          '& ZMassReco < 116 ' \
          '& (ZDecayMode == 13 | ZDecayMode == 151313) ' \
          '& ZLeptonRecoPt > 27 ' \
          '& ZAntiLeptonRecoPt > 27 ' \
          '& abs(ZLeptonRecoEta) < 2.4 ' \
          '& abs(ZAntiLeptonRecoEta) < 2.4 '

ZReco = [tree_to_df(root2array(i, treeName[0], selection=selection, branches=branches),5) for i in inputs]
ZReco = pd.concat(ZReco)

ZReco['eventWeight'] = ZReco['eventWeight']/abs(ZReco['eventWeight'])
nReco = np.sum(ZReco['eventWeight'])

print("inclusive Z Efficiency is: eff = "+str(nReco)+"/"+str(nAcc)+" = "+str(nReco/nAcc))

#diffrential efficiencies

dfOut = pd.DataFrame()


def differentialEff(bins, obs, cuts, sel_reco='', sel_gen=''):
    if sel_reco != '':
        sel_reco = '&'+sel_reco
    if sel_gen != '':
        sel_gen = '&' + sel_gen

    ZEff = [np.zeros(len(bins)-1), np.zeros(len(bins)-1)]

    for i in range(0, len(bins)-1):
        bin_low = bins[i]
        bin_high = bins[i+1]

        nReco_i = np.sum(ZReco.query('{0}>{1} & {0}<{2} {3}'.format(obs, bin_low, bin_high, sel_reco))['eventWeight'])
        nAcc_i  = np.sum(ZAcc.query('{0}>{1} & {0}<{2} {3}'.format(obs, bin_low, bin_high, sel_gen))['eventWeight'])
        if nAcc_i < 100 or nReco < 50:
            continue    # skip low statistic
        ZEff[0][i] = nReco_i/nAcc_i
        ZEff[1][i] = ZEff[0][i] * np.sqrt(1/nReco_i + 1/nAcc_i)

    dfOut[obs + '_' + cuts + '_lowEdge'] = bins[:-1]
    dfOut[obs + '_' + cuts + '_upEdge'] = bins[1:]
    dfOut['ZEff_{0}_{1}'.format(obs, cuts)] = ZEff[0]
    dfOut['ZErr_{0}_{1}'.format(obs, cuts)] = ZEff[1]

    plt.clf()
    fig, ax = plt.subplots()
    ax.errorbar(bins[:-1] + (bins[1:] - bins[:-1])/2, ZEff[0], xerr=np.zeros(len(bins)-1), yerr=ZEff[1], fmt='bo')
    ax.set(xlim=(bins[0], bins[-1]), ylim=(0.7, 1.1))
    ax.text(0.05, 0.9, cuts, transform=ax.transAxes)
    ax.set_ylabel('Z efficiency')
    ax.set_xlabel(obs)
    plt.savefig('ZEff_{0}_{1}.png'.format(obs, cuts))
    plt.close()


differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'inclusive')
differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'BB',
                sel_reco='abs(ZLeptonRecoEta) < 0.9 & abs(ZAntiLeptonRecoEta) < 0.9',
                sel_gen='abs(ZLeptonEta) < 0.9 & abs(ZAntiLeptonEta) < 0.9')
differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'EB',
                sel_reco='abs(ZLeptonRecoEta) > 0.9 & abs(ZAntiLeptonRecoEta) < 0.9',
                sel_gen='abs(ZLeptonEta) > 0.9 & abs(ZAntiLeptonEta) < 0.9')
differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'BE',
                sel_reco='abs(ZLeptonRecoEta) < 0.9 & abs(ZAntiLeptonRecoEta) > 0.9',
                sel_gen='abs(ZLeptonEta) < 0.9 & abs(ZAntiLeptonEta) > 0.9')
differentialEff(np.linspace(0.5, 59.5, 60), 'nPV', 'EE',
                sel_reco='abs(ZLeptonRecoEta) > 0.9 & abs(ZAntiLeptonRecoEta) > 0.9',
                sel_gen='abs(ZLeptonEta) > 0.9 & abs(ZAntiLeptonEta) > 0.9')

dfOut.to_hdf('ZEfficiencies.h5', 'ZEfficiencies', mode='w')



