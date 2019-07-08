from __future__ import division, print_function
# # #
#
#  Measures the true Z->mumu reconstruction efficiency from MC
#
# # #
from pandas import read_hdf
import pdb
import argparse
import numpy as np
import matplotlib.pyplot as plt

dfZEff = read_hdf('ZEfficiencies.h5')
dfZMuMuEff = read_hdf('ZMuMuEfficiencies.h5')

minPU=min(dfZEff['nPV'])
maxPU=max(dfZEff['nPV'])+1

plt.errorbar(range(minPU,maxPU), dfZEff['ZEff'], xerr=np.zeros(maxPU-minPU), yerr=dfZEff['ZErr'], fmt='bo', label='true' )
plt.errorbar(range(minPU,maxPU), dfZMuMuEff['ZMuMuEff_old'], xerr=np.zeros(maxPU-minPU), yerr=dfZMuMuEff['ZMuMuErr_old'], fmt='ro', label='old' )
plt.errorbar(range(minPU,maxPU), dfZMuMuEff['ZMuMuEff_new'], xerr=np.zeros(maxPU-minPU), yerr=dfZMuMuEff['ZMuMuErr_new'], fmt='go', label='new' )
plt.ylim((0.7,1.1))
plt.legend()
plt.ylabel('Z efficiency')
plt.xlabel('nPV')
plt.savefig('ZMuMuEff_PU_all.pdf')
plt.cla()

plt.errorbar(range(minPU,maxPU), dfZMuMuEff['ZMuMuEff_old'] - dfZEff['ZEff'], xerr=np.zeros(maxPU-minPU), yerr=np.sqrt(dfZMuMuEff['ZMuMuErr_old']**2 + dfZEff['ZErr']**2), fmt='ro', label='old' )
plt.errorbar(range(minPU,maxPU), dfZMuMuEff['ZMuMuEff_new'] - dfZEff['ZEff'], xerr=np.zeros(maxPU-minPU), yerr=np.sqrt(dfZMuMuEff['ZMuMuErr_new']**2 + dfZEff['ZErr']**2), fmt='go', label='new' )
plt.ylim((0.,0.2))
plt.legend(loc=2)
plt.ylabel('Z efficiency pulls')
plt.xlabel('nPV')
plt.savefig('ZMuMuEff_PU_pulls.pdf')