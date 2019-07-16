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
from scipy.optimize import curve_fit

dfZEff = read_hdf('ZEfficiencies.h5')
dfZMuMuEff = read_hdf('ZMuMuEfficiencies.h5')


def quadratic(x, a, b, c):
    return a * x**2 + b * x + c


def exponential(x, a, b, c):
    return a * np.exp(x * b) + c


def plotComparisons(obs, cuts):

    bins_low = dfZEff[obs+'_'+cuts+'_lowEdge'].values
    bins_high = dfZEff[obs+'_'+cuts+'_upEdge'].values

    x = bins_low + (bins_high - bins_low)/2

    plt.clf()
    fig, ax = plt.subplots()

    ax.errorbar(x, dfZEff['ZEff_{0}_{1}'.format(obs, cuts)], xerr=np.zeros(len(x)), yerr=dfZEff['ZErr_{0}_{1}'.format(obs, cuts)], fmt='bo', label='true' )
    ax.errorbar(x, dfZMuMuEff['ZMuMuEff_{0}_{1}_old'.format(obs, cuts)], xerr=np.zeros(len(x)), yerr=dfZMuMuEff['ZMuMuErr_{0}_{1}_old'.format(obs, cuts)], fmt='ro', label='old' )
    ax.errorbar(x, dfZMuMuEff['ZMuMuEff_{0}_{1}_new'.format(obs, cuts)], xerr=np.zeros(len(x)), yerr=dfZMuMuEff['ZMuMuErr_{0}_{1}_new'.format(obs, cuts)], fmt='go', label='new' )
    ax.set(xlim=(bins_low[0], bins_high[-1]), ylim=(0.7, 1.1))
    ax.legend()
    ax.text(0.05, 0.9, cuts, transform=ax.transAxes)
    ax.set_ylabel('Z efficiency')
    ax.set_xlabel(obs)
    plt.savefig('ZMuMu_{0}_{1}_all.png'.format(obs, cuts))
    plt.close()



    pulls_old = dfZMuMuEff['ZMuMuEff_{0}_{1}_old'.format(obs, cuts)] - dfZEff['ZEff_{0}_{1}'.format(obs, cuts)]
    pulls_old_sig = np.sqrt(dfZMuMuEff['ZMuMuErr_{0}_{1}_old'.format(obs, cuts)]**2 + dfZEff['ZErr_{0}_{1}'.format(obs, cuts)]**2)
    pulls_new = dfZMuMuEff['ZMuMuEff_{0}_{1}_new'.format(obs, cuts)] - dfZEff['ZEff_{0}_{1}'.format(obs, cuts)]
    pulls_new_sig = np.sqrt(dfZMuMuEff['ZMuMuErr_{0}_{1}_new'.format(obs, cuts)]**2 + dfZEff['ZErr_{0}_{1}'.format(obs, cuts)]**2)

    start = 5
    stop = 45

    popt, pcov = curve_fit(quadratic, x[start:stop], pulls_old[start:stop], sigma=pulls_old_sig[start:stop], absolute_sigma=True)
    popt_exp, pcov_exp = curve_fit(exponential, x[start:stop], pulls_old[start:stop], p0=np.asarray([0.001, 0.001, 0.01]), sigma=pulls_old_sig[start:stop], absolute_sigma=True)

    plt.clf()
    fig, ax = plt.subplots()
    ax.errorbar(x, pulls_old, xerr=np.zeros(len(x)), yerr=pulls_old_sig, fmt='ro', label='old')
    ax.errorbar(x, pulls_new, xerr=np.zeros(len(x)), yerr=pulls_new_sig, fmt='go', label='new')
    ax.plot(x, quadratic(x, *popt), 'c--', label='{0:.2e} x^2 + {1:.2e} x + {2:.2e}'.format(popt[0], popt[1], popt[2]), linewidth=2)
    ax.plot(x, exponential(x, *popt_exp), 'm:', label='{0:.2e} * exp({1:.2e} * x) + {2:.2e}'.format(popt_exp[0], popt_exp[1], popt_exp[2]), linewidth=2)
    ax.set(xlim=(bins_low[0], bins_high[-1]), ylim=(0.0, 0.1))
    ax.text(0.05, 0.625, cuts, transform=ax.transAxes)
    ax.legend(loc=2)
    ax.set_ylabel('Z efficiency pulls')
    ax.set_xlabel(obs)
    plt.savefig('ZMuMu_{0}_{1}_pulls.png'.format(obs, cuts))
    plt.close()


plotComparisons('nPV', 'inclusive')
plotComparisons('nPV', 'BB')
plotComparisons('nPV', 'BE')
plotComparisons('nPV', 'EB')
plotComparisons('nPV', 'EE')
