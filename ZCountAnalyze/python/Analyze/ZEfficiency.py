from __future__ import division, print_function
# # #
#
#  Measures the true Z->mumu reconstruction efficiency from MC
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

parser = argparse.ArgumentParser(prog='./ZEfficiency')
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

#dataframe with all Zs that are accepted in Gen
ZAcc = tree_to_df(root2array(input, treeName[0], selection=selection, branches=branches),5)
ZAcc['eventWeight'] = ZAcc['eventWeight']/abs(ZAcc['eventWeight'])
nAcc = np.sum(ZAcc['eventWeight'])

print("number of accepted gen Events is = "+str(nAcc))


#dataframe with all Zs that are reconstructed
selection='nMuonTags >=1 & nMuonProbes >= 1'
ZReco = tree_to_df(root2array(input, treeName[0], selection=selection, branches=branches),5)
query = ' | '
query = query.join(['MuonProbeCategory_{0}==1 | MuonProbeCategory_{0}==2'.format(i) for i in range(0,5)])
ZReco = ZReco.query(query)
ZReco['eventWeight'] = ZReco['eventWeight']/abs(ZReco['eventWeight'])
nReco = np.sum(ZReco['eventWeight'])

print("inclusive Z Efficiency is: eff = "+str(nReco)+"/"+str(nAcc)+" = "+str(nReco/nAcc))

#Efficiency in dependence of PU
minPU=0
maxPU=60

ZEff_PU = np.zeros(maxPU-minPU)
ZEff_PU_err = np.zeros(maxPU-minPU)
for iPV in range(minPU,maxPU):
    nReco_i = np.sum(ZReco.query('nPV=={0}'.format(iPV))['eventWeight'])
    nAcc_i  = np.sum(ZAcc.query('nPV=={0}'.format(iPV))['eventWeight'])
    if nAcc_i < 100 or nReco < 50:
        continue #skip low statistic
    ZEff_PU[iPV] = nReco_i/nAcc_i
    ZEff_PU_err[iPV] = ZEff_PU[iPV] *np.sqrt( 1/nReco_i + 1/nAcc_i )


dfOut = DataFrame()
dfOut['nPV'] = range(minPU,maxPU)
dfOut['ZEff'] = ZEff_PU
dfOut['ZErr'] = ZEff_PU_err
dfOut.to_hdf('ZEfficiencies.h5', 'ZEfficiencies', mode='w')


### Plot ###

plt.errorbar(range(minPU,maxPU), ZEff_PU, xerr=np.zeros(maxPU-minPU), yerr=ZEff_PU_err, fmt='bo' )
plt.ylim((0.7,1.1))
plt.ylabel('Z efficiency')
plt.xlabel('nPV')
plt.savefig('ZEff_PU.pdf')

