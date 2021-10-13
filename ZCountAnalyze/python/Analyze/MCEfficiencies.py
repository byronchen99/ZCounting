from __future__ import division, print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
from root_numpy import root2array, list_trees
import math
import json
import pdb
import uncertainties as unc
import argparse

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df

massLo=56
massHi=116
ptCut=30

def zMCEff1D(dfGen, dfReco, bins, obs, region='inclusive', selGen='', selReco='',
             cutsPtEta='p_\mathrm{T}(\mu) > '+str(ptCut)+'\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4'):
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

    if selGen != '':
        dfGen = dfGen.query(selGen)
    if selReco != '':
        dfReco = dfReco.query(selReco)

    Efficiency.Set(bins, obs, region)

    for i in range(0, len(bins) - 1):
        bin_low = bins[i]
        bin_high = bins[i + 1]
        sel = '{0}>{1} & {0}<{2}'.format(obs, bin_low, bin_high)

        Efficiency.Calculate(dfGen.query(sel), dfReco.query(sel), i)

    Efficiency.Plot(cutsPtEta)


class Efficiency:
    collection = []

    @classmethod
    def Set(cls, bins, obs, region):
        for eff in cls.collection:
            eff.set(bins, obs, region)

    @classmethod
    def Calculate(cls, dfGen, dfReco, ibin):
        for eff in cls.collection:
            eff.calculate(dfGen, dfReco, ibin)

    @classmethod
    def Plot(cls, cutsPtEta):
        for eff in cls.collection:
            eff.plot(cutsPtEta)

    def __init__(self, name, reqID, reqPFIso=None, reqTkIso=None,
        delR=None, delDz=None, delDxy=None,
        reqHLT=None, triggerName='', range=(0.9,1.0)):
        self.name = name
        self.ID = reqID
        self.pfIso = reqPFIso
        self.tkIso = reqTkIso
        self.delR = delR
        self.delDz = delDz
        self.delDxy = delDxy
        self.HLT = reqHLT
        self.triggerName = triggerName
        self.eff_corr = []
        self.eff_true = []
        self.eff_tnp = []
        self.eff_tnpZ = []
        self.bins = []
        self.obs = None
        self.region = None
        self.x = []
        self.range = range

        self.fitResults = dict()
        Efficiency.collection.append(self)

    def set(self, bins, obs, region):
        self.bins = bins
        self.obs = obs
        self.region = region
        self.x = bins[:-1] + (bins[1:] - bins[:-1]) / 2
        self.eff_corr = []
        self.eff_true = []
        self.eff_tnp = []
        self.eff_tnpZ = []

    def calculate(self, dataGen, dataReco, ibin):

        # >>> acceptance selection
        # var0 = [True,]*len(dataframeReco)
        # if self.delR:
        #     var0 = var0 & (dataframe['delR'] > self.delR)
        # if self.delDxy:
        #     var0 = var0 & (dataframe['delDxy'] < self.delDxy)
        # if self.delDz:
        #     var0 = var0 & (dataframe['delDz'] < self.delDz)

        # if all(var0):
        #     dfGen = dataframeReco.copy()
        # else:
        #     dfGen = dataframeReco.loc[var0].copy()
        # <<< acceptance selection

        # >>> columns for tnp efficiency
        dfReco = dataReco.copy()
        dfReco['pass_muon'] = dfReco['muon_ID'] >= self.ID
        dfReco['pass_antiMuon'] = dfReco['antiMuon_ID'] >= self.ID
        if self.pfIso:
            dfReco['pass_muon'] = dfReco['pass_muon'] & (dfReco['muon_pfIso'] < self.pfIso)
            dfReco['pass_antiMuon'] = dfReco['pass_antiMuon'] & (dfReco['antiMuon_pfIso'] < self.pfIso)

        if self.tkIso:
            dfReco['pass_muon'] = dfReco['pass_muon'] & (dfReco['muon_tkIso'] < self.tkIso)
            dfReco['pass_antiMuon'] = dfReco['pass_antiMuon'] & (dfReco['antiMuon_tkIso'] < self.tkIso)

        dfReco['pass_muon_antiMuon'] = dfReco['pass_muon'] & dfReco['pass_antiMuon']

        nMinus = dfReco['pass_muon'].sum()
        nPlus = dfReco['pass_antiMuon'].sum()
        nPlusMinus = dfReco['pass_muon_antiMuon'].sum()

        if self.HLT:
            dfReco['pass_muon_HLT'] = dfReco['pass_muon_antiMuon'] & dfReco['muon_hlt_{0}'.format(self.HLT)]
            dfReco['pass_antiMuon_HLT'] = dfReco['pass_muon_antiMuon'] & dfReco['antiMuon_hlt_{0}'.format(self.HLT)]
            dfReco['pass_muon_antiMuon_HLT'] = dfReco['pass_muon_HLT'] & dfReco['pass_antiMuon_HLT']

            nMinus_HLT = dfReco['pass_muon_HLT'].sum()
            nPlus_HLT = dfReco['pass_antiMuon_HLT'].sum()
            nPlusMinus_HLT = dfReco['pass_muon_antiMuon_HLT'].sum()
        # <<< columns for tnp efficiency

        # >>> columns for true efficiency
        dfGen = dataGen.copy()
        dfGen['pass_muon'] = dfGen['muon_ID'] >= self.ID
        dfGen['pass_antiMuon'] = dfGen['antiMuon_ID'] >= self.ID
        if self.pfIso:
            dfGen['pass_muon'] = dfGen['pass_muon'] & (dfGen['muon_pfIso'] < self.pfIso)
            dfGen['pass_antiMuon'] = dfGen['pass_antiMuon'] & (dfGen['antiMuon_pfIso'] < self.pfIso)

        if self.tkIso:
            dfGen['pass_muon'] = dfGen['pass_muon'] & (dfGen['muon_tkIso'] < self.tkIso)
            dfGen['pass_antiMuon'] = dfGen['pass_antiMuon'] & (dfGen['antiMuon_tkIso'] < self.tkIso)

        dfGen['pass_muon_antiMuon'] = dfGen['pass_muon'] & dfGen['pass_antiMuon']

        if self.HLT:
            dfGen['pass_muon_HLT'] = dfGen['pass_muon_antiMuon'] & dfGen['muon_hlt_{0}'.format(self.HLT)]
            dfGen['pass_antiMuon_HLT'] = dfGen['pass_muon_antiMuon'] & dfGen['antiMuon_hlt_{0}'.format(self.HLT)]

            dfGen['reco'] = dfGen['pass_muon_antiMuon'] & (dfGen['pass_muon_HLT'] | dfGen['pass_antiMuon_HLT'])
        else:
            dfGen['reco'] = dfGen['pass_muon_antiMuon']

        dfGen['not_reco'] = ~dfGen['reco']

        nZReco = dfGen['reco'].sum()
        nZNotReco = dfGen['not_reco'].sum()
        # <<< columns for true efficiency

        # # >>> construct full covariance matrix
        # if self.HLT:
        #     corr_matrix = dfReco[['reco','not_reco','pass_muon_antiMuon','pass_muon','pass_antiMuon','pass_muon_antiMuon_HLT','pass_muon_HLT','pass_antiMuon_HLT']].corr()
        #
        #     nZReco, nZNotReco, nPlusMinus, nMinus, nPlus, nPlusMinus_HLT, nMinus_HLT, nPlus_HLT = unc.correlated_values_norm(
        #         [(nZReco, np.sqrt(nZReco)), (nZNotReco, np.sqrt(nZNotReco)),
        #         (nPlusMinus, np.sqrt(nPlusMinus)), (nMinus, np.sqrt(nMinus)), (nPlus, np.sqrt(nPlus)),
        #         (nPlusMinus_HLT, np.sqrt(nPlusMinus_HLT)), (nMinus_HLT, np.sqrt(nMinus_HLT)), (nPlus_HLT, np.sqrt(nPlus_HLT))],
        #         corr_matrix)
        #
        # else:
        #     corr_matrix = dfReco[['reco','not_reco','pass_muon_antiMuon','pass_muon','pass_antiMuon']].corr()
        #
        #     nZReco, nZNotReco, nPlusMinus, nMinus, nPlus = unc.correlated_values_norm(
        #         [(nZReco, np.sqrt(nZReco)), (nZNotReco, np.sqrt(nZNotReco)),
        #         (nPlusMinus, np.sqrt(nPlusMinus)), (nMinus, np.sqrt(nMinus)), (nPlus, np.sqrt(nPlus))],
        #         corr_matrix)
        # # <<< construct full covariance matrix

        # >>> compute true and tnp efficiency
        eff_true = nZReco / (nZNotReco + nZReco)

        MuPlusEff = nPlusMinus/nMinus
        MuMinusEff = nPlusMinus/nPlus

        eff_tnp = MuPlusEff * MuMinusEff

        if self.HLT:
            MuPlusEff_HLT = nPlusMinus_HLT/nMinus_HLT
            MuMinusEff_HLT = nPlusMinus_HLT/nPlus_HLT

            eff_tnpZ = (1 - (1 - MuPlusEff_HLT) * (1 - MuMinusEff_HLT) ) * eff_tnp
        # <<< compute true and tnp efficiency
        else:
            eff_tnpZ = eff_tnp

        # >>> store
        # self.eff_corr.append(unc.covariance_matrix([eff_tnpZ, eff_true]))
        self.eff_tnp.append(eff_tnp)
        self.eff_tnpZ.append(eff_tnpZ)
        self.eff_true.append(eff_true)



    def plot(self, cutsPtEta):
        eff_tnp = self.eff_tnp # np.array([x.nominal_value for x in self.eff_tnp])
        eff_tnpZ = self.eff_tnpZ # np.array([x.nominal_value for x in self.eff_tnpZ])
        # eff_tnpZ_err = np.array([x.std_dev for x in self.eff_tnpZ])
        eff_true = self.eff_true # np.array([x.nominal_value for x in self.eff_true])
        # eff_true_err = np.array([x.std_dev for x in self.eff_true])

        plt.clf()
        fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        fig.subplots_adjust(hspace=0)
        fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.1)
        ax[0].plot(self.x, eff_tnp,  color='grey', linestyle='-', label=r'$\mathrm{TnP}\ (\epsilon^\mu_\mathrm{Reco})^2$')
        ax[0].plot(self.x, eff_tnpZ,  color='blue', linestyle='-', label=r'$\mathrm{TnP}\ \epsilon^\mathrm{Z}_\mathrm{Reco}$')
        ax[0].plot(self.x, eff_true, color='red', linestyle='--', label=r'$\mathrm{True}\ \epsilon^\mathrm{Z}_\mathrm{Reco}$')

        # ax[0].errorbar(self.x, eff_true, xerr=np.zeros(len(self.x)), yerr=np.zeros(len(eff_true_err)), #eff_true_err,
        #     fmt='-b', label='true', markeredgecolor='black')
        # ax[0].errorbar(self.x, eff_tnpZ, xerr=np.zeros(len(self.x)), yerr=np.zeros(len(eff_tnpZ_err)), #eff_tnpZ_err,
        #     fmt='--r', label='tnp', markeredgecolor='black')
        ymin = min([min(eff_true), min(eff_tnpZ)])
        ymax = max([max(eff_true), max(eff_tnpZ)])
        ymax += (ymax-ymin)*0.4
        ymin -= (ymax-ymin)*0.05
        # ax[0].set(xlim=(self.bins[0], self.bins[-1]), ylim=(ymin, ymax))
        ax[0].set(xlim=(self.bins[0], self.bins[-1]), ylim=self.range)
        ax[0].legend(frameon=False)
        # legend = ax[0].legend()
        # legend.get_frame().set_facecolor('none')
        if self.delR:
            ax[0].text(0.5, 0.9, r'$\Delta R(\mu,\mu) > '+str(self.delR)+"$",
                        transform=ax[0].transAxes)
        if self.delDxy:
            ax[0].text(0.5, 0.82, r'$\Delta D_\mathrm{XY}(\mu,\mu) < '+str(self.delDxy)+"$",
                        transform=ax[0].transAxes)
        if self.delDz:
            ax[0].text(0.5, 0.74, r'$\Delta D_\mathrm{Z}(\mu,\mu) < '+str(self.delDz)+"$",
                        transform=ax[0].transAxes)

        if self.pfIso:
            ax[0].text(0.05, 0.28, r'$Iso_\mathrm{PF}(\mu) < '+str(self.pfIso)+"$",
                        transform=ax[0].transAxes)
        if self.tkIso:
            ax[0].text(0.05, 0.20, r'$Iso_\mathrm{Tk}(\mu) < '+str(self.tkIso)+"$",
                        transform=ax[0].transAxes)

        if self.ID == 4:
            ax[0].text(0.05, 0.12, r'tight ID (w/o IPC)', transform=ax[0].transAxes)
        elif self.ID == 5:
            ax[0].text(0.05, 0.12, r'tight ID', transform=ax[0].transAxes)

        ax[0].text(0.05, 0.9, str(massLo)+r'$\mathrm{GeV} < m(\mu,\mu) < '+str(massHi)+r'\ \mathrm{GeV}$',
            transform=ax[0].transAxes)
        ax[0].text(0.05, 0.82, r'${0}$'.format(cutsPtEta),
            transform=ax[0].transAxes)

        if self.HLT:
            ax[0].text(0.05, 0.04,
                       r'{0}'.format(self.triggerName),
                       transform=ax[0].transAxes, color='g')
        ax[0].set_ylabel(r'$\epsilon_\mathrm{Z}$')
        ax[0].set_xlabel(self.obs)
        # ax[0].set_yticks([0.8, 0.85, 0.9, 0.95, 1.0, 1.05])
        #ax[0].set_title("Z efficiency at {0} muon level".format(self.name))

        pulls = [x - y for (x,y) in zip(self.eff_tnpZ, self.eff_true)]
        pulls_val = pulls#[x.nominal_value for x in pulls]
        #pulls_err = [x.std_dev for x in pulls]

        # do fit
        def linear(x, a, b):
            return a * x + b

        popt, pcov = curve_fit(linear, self.x, pulls_val)#, sigma=pulls_err)
        perr = np.sqrt(np.diag(pcov))
        self.fitResults[self.region + "_a"] = popt[0]
        self.fitResults[self.region + "_b"] = popt[1]
        self.fitResults[self.region + "_a-err"] = perr[0]
        self.fitResults[self.region + "_b-err"] = perr[1]

        ax[1].plot(self.x, np.zeros(len(self.x)), color='gray', linestyle='dashed')

        ax[1].plot(self.x, linear(self.x, *popt), 'm-')
        ax[1].fill_between(self.x, linear(self.x, *(popt-perr)), linear(self.x, *(popt+perr)), color='magenta', alpha=0.5)
        if popt[1] > 0:
            ax[1].text(0.04, 0.05, r'fit: $(%5.1e\ \pm\ %5.1e) \cdot x + %5.1e\ \pm\ %5.1e $' % (popt[0],perr[0], popt[1],perr[1]), color='black', transform=ax[1].transAxes)
        else:
            ax[1].text(0.04, 0.05, r'fit: $(%5.1e\ \pm\ %5.1e) \cdot x - %5.1e\ \pm\ %5.1e $' % (popt[0],perr[0], abs(popt[1]),perr[1]), color='black', transform=ax[1].transAxes)

        ax[1].errorbar(self.x, pulls_val, xerr=np.zeros(len(self.x)), # yerr=pulls_err,
                       fmt='ko')  # , label='factorized - true')
        ax[1].set_ylim(-0.03, 0.03)
        ax[1].set_ylabel(r'$\epsilon^\mathrm{tnp}_\mathrm{Z} - \epsilon^\mathrm{true}_\mathrm{Z}$')
        if self.obs == 'nPU':
            ax[1].set_xlabel(r'$N^\mathrm{PU}$')
        else:
            ax[1].set_xlabel(self.obs)
        ax[1].set_yticks([-0.02, -0.01, 0., 0.01, 0.02])
        # fig.tight_layout()
        plt.savefig(output + '/ZMuMu_{0}_{1}_{2}_level.eps'.format(self.obs, self.region, self.name))
        plt.savefig(output + '/ZMuMu_{0}_{1}_{2}_level.png'.format(self.obs, self.region, self.name))
        plt.close()


# acceptance selection
selectionGen = 'decayMode == 13 ' \
            '& muon_genRecoMatches == 1' \
            '& antiMuon_genRecoMatches == 1' \

#             '& z_genMass > {0} ' \
#             '& z_genMass < {1} ' \
#             '& muon_genPt > {2} ' \
#             '& antiMuon_genPt > {2} ' \
#             '& abs(muon_genEta) < 2.4 ' \
#             '& abs(antiMuon_genEta) < 2.4 ' \
#             '& muon_genRecoMatches == 1' \
#             '& antiMuon_genRecoMatches == 1'.format(massLo, massHi, ptCut)

selectionReco = 'z_recoMass > {0} ' \
            '& z_recoMass < {1} ' \
            '& muon_genRecoMatches == 1' \
            '& antiMuon_genRecoMatches == 1' \
            '& decayMode == 13'.format(massLo, massHi)

# specify which branches to load
branches = [
            # 'nPV', 'nPU',
            'eventweight',
            # 'z_genMass',
            # 'z_recoMass',
            # 'muon_genEta', 'antiMuon_genEta',
            # 'muon_genPhi', 'antiMuon_genPhi',
            # 'muon_genPt', 'antiMuon_genPt',
            # 'muon_genRecoMatches',  'antiMuon_genRecoMatches',
            # 'muon_genRecoObj',      'antiMuon_genRecoObj',
            # 'Muon_dxy[muon_genRecoObj]',            'Muon_dxy[antiMuon_genRecoObj]',
            # 'Muon_dz[muon_genRecoObj]',             'Muon_dz[antiMuon_genRecoObj]',
            'Muon_pt[muon_genRecoObj]',             'Muon_pt[antiMuon_genRecoObj]',
            'Muon_eta[muon_genRecoObj]',            'Muon_eta[antiMuon_genRecoObj]',
            'Muon_phi[muon_genRecoObj]',            'Muon_phi[antiMuon_genRecoObj]',
            # 'Muon_ID[muon_genRecoObj]',             'Muon_ID[antiMuon_genRecoObj]',
            # 'Muon_triggerBits[muon_genRecoObj]',    'Muon_triggerBits[antiMuon_genRecoObj]',
            # 'Muon_tkRelIso[muon_genRecoObj]',       'Muon_tkRelIso[antiMuon_genRecoObj]',
            # 'Muon_pfRelIso04_all[muon_genRecoObj]', 'Muon_pfRelIso04_all[antiMuon_genRecoObj]',
            ]


parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-i', '--input', nargs='+',
    help='specify input root file'
)
parser.add_argument(
    '-y', '--year', default=2017, type=int,
    help='specify year'
)
parser.add_argument(
    '-o', '--output', nargs=1, default='./',
    help='specify output dir'
)
args = parser.parse_args()

inputs = args.input
output = args.output[0]

if os.path.isfile(output+"/dfGen.h5") and os.path.isfile(output+"/dfReco.h5"):
    print(">>> load dataframe")
    store = pd.HDFStore(output+'/dfGen.h5')
    dfGen = store['dfGen']  # load it
    store = pd.HDFStore(output+'/dfReco.h5')
    dfReco = store['dfReco']  # load it
else:
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

    print(">>> Load Events in gen acceptance")
    # >>> in fiducial phase space of Gen info for true efficiency
    dfGen = [tree_to_df(root2array(i, treeName[0], selection=selectionGen, branches=branches), 1) for i in inputs]
    dfGen = pd.concat(dfGen)

    dfGen = dfGen.rename(columns={
        'Muon_ID[muon_genRecoObj]_0': 'muon_ID',
        'Muon_ID[antiMuon_genRecoObj]_0': 'antiMuon_ID',
        'Muon_triggerBits[muon_genRecoObj]_0': 'muon_triggerBits',
        'Muon_triggerBits[antiMuon_genRecoObj]_0': 'antiMuon_triggerBits',
        # 'Muon_tkRelIso[muon_genRecoObj]_0': 'muon_tkIso',
        # 'Muon_tkRelIso[antiMuon_genRecoObj]_0': 'antiMuon_tkIso',
        # 'Muon_pfRelIso04_all[muon_genRecoObj]_0': 'muon_pfIso',
        # 'Muon_pfRelIso04_all[antiMuon_genRecoObj]_0': 'antiMuon_pfIso',
        'Muon_eta[muon_genRecoObj]_0': 'muon_eta',
        'Muon_eta[antiMuon_genRecoObj]_0': 'antiMuon_eta',
        'Muon_pt[muon_genRecoObj]_0': 'muon_pt',
        'Muon_pt[antiMuon_genRecoObj]_0': 'antiMuon_pt',
        'Muon_phi[muon_genRecoObj]_0': 'muon_phi',
        'Muon_phi[antiMuon_genRecoObj]_0': 'antiMuon_phi',
        # 'Muon_dxy[muon_genRecoObj]_0': 'muon_dxy',
        # 'Muon_dxy[antiMuon_genRecoObj]_0': 'antiMuon_dxy',
        # 'Muon_dz[muon_genRecoObj]_0': 'muon_dz',
        # 'Muon_dz[antiMuon_genRecoObj]_0': 'antiMuon_dz'
        })

    print(">>> store dataframe")
    store = pd.HDFStore(output+'/dfGen.h5')
    store['dfGen'] = dfGen  # save it
    exit()

    # print(">>> add new columns")
    #dfGen['delPhiLL'] = abs(
    #    abs(dfGen['muon_genPhi'] - dfGen['antiMuon_genPhi']).apply(lambda x: x - 2 * math.pi if x > math.pi else x))
    #dfGen['delEtaLL'] = abs(dfGen['muon_genEta'] - dfGen['antiMuon_genEta'])
    # dfGen['delR'] = np.sqrt(
    #    (dfGen['muon_eta'] - dfGen['antiMuon_eta']) ** 2 + (dfGen['muon_phi'] - dfGen['antiMuon_phi']) ** 2)
    #
    # dfGen['delDz'] = abs(dfGen['muon_dz'] - dfGen['antiMuon_dz'])
    # dfGen['delDxy'] = abs(dfGen['muon_dxy'] - dfGen['antiMuon_dxy'])

    #dfGen['delPtLL'] = abs(dfGen['muon_genPt'] - dfGen['antiMuon_genPt'])
    #dfGen['relPtLL'] = abs(dfGen['muon_genPt'] - dfGen['antiMuon_genPt']) / abs(
    #    dfGen['muon_genPt'] + dfGen['antiMuon_genPt'])
    #dfGen['sumPtLL'] = dfGen['muon_genPt'] + dfGen['antiMuon_genPt']

    #dfGen = dfGen.query('delRLL > 0.4')

    print(">>> convert bit code into bit map")
    #   1: "HLT_L1SingleMu18_v*"
    #   2: "HLT_L1SingleMu25_v*"
    #   3: "HLT_IsoMu24_v*"
    #   4: "HLT_IsoTkMu24_v*"
    #   5: "HLT_IsoMu27_v*"
    #   6: "HLT_IsoTkMu27_v*"
    #   7: "HLT_IsoMu30_v*"
    #   8: "HLT_IsoTkMu30_v*"

    for iBit in range(0, 10):
        nBit = 2 ** iBit
        dfGen['muon_hlt_{0}'.format(iBit+1)] = dfGen['muon_triggerBits'].apply(
            lambda x: 1 if x % (nBit * 2) >= nBit else 0)
        dfGen['antiMuon_hlt_{0}'.format(iBit+1)] = dfGen['antiMuon_triggerBits'].apply(
            lambda x: 1 if x % (nBit * 2) >= nBit else 0)


    print(">>> store dataframe")
    store = pd.HDFStore(output+'/dfGen.h5')
    store['dfGen'] = dfGen  # save it

    # >>> in fiducial phase space of Reco info for true efficiency
    print(">>> Load Events in reco acceptance")
    dfReco = [tree_to_df(root2array(i, treeName[0], selection=selectionReco, branches=branches), 1) for i in inputs]
    dfReco = pd.concat(dfReco)

    dfReco = dfReco.rename(columns={
        'Muon_ID[muon_genRecoObj]_0': 'muon_ID',
        'Muon_ID[antiMuon_genRecoObj]_0': 'antiMuon_ID',
        'Muon_triggerBits[muon_genRecoObj]_0': 'muon_triggerBits',
        'Muon_triggerBits[antiMuon_genRecoObj]_0': 'antiMuon_triggerBits',
        # 'Muon_tkRelIso[muon_genRecoObj]_0': 'muon_tkIso',
        # 'Muon_tkRelIso[antiMuon_genRecoObj]_0': 'antiMuon_tkIso',
        # 'Muon_pfRelIso04_all[muon_genRecoObj]_0': 'muon_pfIso',
        # 'Muon_pfRelIso04_all[antiMuon_genRecoObj]_0': 'antiMuon_pfIso',
        'Muon_eta[muon_genRecoObj]_0': 'muon_eta',
        'Muon_eta[antiMuon_genRecoObj]_0': 'antiMuon_eta',
        'Muon_pt[muon_genRecoObj]_0': 'muon_pt',
        'Muon_pt[antiMuon_genRecoObj]_0': 'antiMuon_pt',
        # 'Muon_phi[muon_genRecoObj]_0': 'muon_phi',
        # 'Muon_phi[antiMuon_genRecoObj]_0': 'antiMuon_phi',
        # 'Muon_dxy[muon_genRecoObj]_0': 'muon_dxy',
        # 'Muon_dxy[antiMuon_genRecoObj]_0': 'antiMuon_dxy',
        # 'Muon_dz[muon_genRecoObj]_0': 'muon_dz',
        # 'Muon_dz[antiMuon_genRecoObj]_0': 'antiMuon_dz'
        })

    dfReco = dfReco.query('muon_pt > {0} & antiMuon_pt > {0} & abs(muon_eta) < 2.4 & abs(antiMuon_eta) < 2.4 '.format(ptCut))

    # print(">>> add new columns")
    #dfReco['delPhiLL'] = abs(
    #    abs(dfReco['muon_genPhi'] - dfReco['antiMuon_genPhi']).apply(lambda x: x - 2 * math.pi if x > math.pi else x))
    #dfReco['delEtaLL'] = abs(dfReco['muon_genEta'] - dfReco['antiMuon_genEta'])
    # dfReco['delR'] = np.sqrt(
    #    (dfReco['muon_eta'] - dfReco['antiMuon_eta']) ** 2 + (dfReco['muon_phi'] - dfReco['antiMuon_phi']) ** 2)
    #
    # dfReco['delDz'] = abs(dfReco['muon_dz'] - dfReco['antiMuon_dz'])
    # dfReco['delDxy'] = abs(dfReco['muon_dxy'] - dfReco['antiMuon_dxy'])

    #dfReco['delPtLL'] = abs(dfReco['muon_genPt'] - dfReco['antiMuon_genPt'])
    #dfReco['relPtLL'] = abs(dfReco['muon_genPt'] - dfReco['antiMuon_genPt']) / abs(
    #    dfReco['muon_genPt'] + dfReco['antiMuon_genPt'])
    #dfReco['sumPtLL'] = dfReco['muon_genPt'] + dfReco['antiMuon_genPt']

    #dfReco = dfReco.query('delRLL > 0.4')

    print(">>> convert bit code into bit map")
    #   1: "HLT_L1SingleMu18_v*"
    #   2: "HLT_L1SingleMu25_v*"
    #   3: "HLT_IsoMu24_v*"
    #   4: "HLT_IsoTkMu24_v*"
    #   5: "HLT_IsoMu27_v*"
    #   6: "HLT_IsoTkMu27_v*"
    #   7: "HLT_IsoMu30_v*"
    #   8: "HLT_IsoTkMu30_v*"

    for iBit in range(0, 10):
        nBit = 2 ** iBit
        dfReco['muon_hlt_{0}'.format(iBit+1)] = dfReco['muon_triggerBits'].apply(
            lambda x: 1 if x % (nBit * 2) >= nBit else 0)
        dfReco['antiMuon_hlt_{0}'.format(iBit+1)] = dfReco['antiMuon_triggerBits'].apply(
            lambda x: 1 if x % (nBit * 2) >= nBit else 0)


    print(">>> store dataframe")
    store = pd.HDFStore(output+'/dfReco.h5')
    store['dfReco'] = dfReco  # save it

# --- define Efficiencies
print(">>> define efficiencies")
if args.year == 2016:
    xLow = 0.5
    xHigh = 74.5
    steps = 15

    dfGen['muon_hlt_1'] = dfGen['muon_hlt_3'] | dfGen['muon_hlt_4']
    dfGen['antiMuon_hlt_1'] = dfGen['antiMuon_hlt_3'] | dfGen['antiMuon_hlt_4']
    dfReco['muon_hlt_1'] = dfReco['muon_hlt_3'] | dfReco['muon_hlt_4']
    dfReco['antiMuon_hlt_1'] = dfReco['antiMuon_hlt_3'] | dfReco['antiMuon_hlt_4']

    eff_ZTightIDIsoMu24OrIsoTkMu24 = Efficiency(
        name="ZcTightID_IsoMu24OrIsoTkMu24",
        reqID=4, reqHLT=1,
        triggerName='$\mathrm{HLT\_Iso(TK)Mu24\_v^{*}}$',
        range=(0.9,1.0))

else:
    xLow = 0.5
    xHigh = 74.5
    steps = 15

    ZcTightID_L1SingleMu25 = Efficiency(
        name="ZcTightID_L1SingleMu25",
        reqID=4, reqHLT=2,
        triggerName='$\mathrm{HLT\_L1SingleMu25\_v^{*}}$'
        )

    ZcTightID_IsoMu27 = Efficiency(
        name="ZcTightID_IsoMu27",
        reqID=4, reqHLT=4,
        triggerName='$\mathrm{HLT\_IsoMu27\_v^{*}}$'
        )


zMCEff1D(dfGen, dfReco, np.linspace(xLow, xHigh, steps), 'nPU', 'inclusive')
zMCEff1D(dfGen, dfReco, np.linspace(xLow, xHigh, steps), 'nPU', 'BB',
         selGen='abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) < 0.9',
         selReco='abs(muon_eta) < 0.9 & abs(antiMuon_eta) < 0.9',
         cutsPtEta='p_\mathrm{T}(\mu) > 30\ \mathrm{GeV} \qquad |\eta(\mu)| < 0.9')

zMCEff1D(dfGen, dfReco, np.linspace(xLow, xHigh, steps), 'nPU', 'BE',
         selGen='(abs(muon_genEta) < 0.9 & abs(antiMuon_genEta) > 0.9) | (abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) < 0.9)',
         selReco='(abs(muon_eta) < 0.9 & abs(antiMuon_eta) > 0.9) | (abs(muon_eta) > 0.9 & abs(antiMuon_eta) < 0.9)')

zMCEff1D(dfGen, dfReco, np.linspace(xLow, xHigh, steps), 'nPU', 'EE',
         selGen='abs(muon_genEta) > 0.9 & abs(antiMuon_genEta) > 0.9',
         selReco='abs(muon_eta) > 0.9 & abs(antiMuon_eta) > 0.9',
         cutsPtEta='p_\mathrm{T}(\mu) > 30\ \mathrm{GeV} \qquad 0.9 < |\eta(\mu)| < 2.4')


if args.year == 2016:
    with open(output+'/MCCorrections_nPU_'+ eff_ZTightIDIsoMu24OrIsoTkMu24.name + '.json', 'w') as file:
        file.write(json.dumps(eff_ZTightIDIsoMu24OrIsoTkMu24.fitResults, sort_keys=True, indent=4))
else:
    with open(output+'/MCCorrections_nPU_'+ ZcTightID_L1SingleMu25.name + '.json', 'w') as file:
        file.write(json.dumps(ZcTightID_L1SingleMu25.fitResults, sort_keys=True, indent=4))

    with open(output+'/MCCorrections_nPU_'+ ZcTightID_IsoMu27.name + '.json', 'w') as file:
        file.write(json.dumps(ZcTightID_IsoMu27.fitResults, sort_keys=True, indent=4))
