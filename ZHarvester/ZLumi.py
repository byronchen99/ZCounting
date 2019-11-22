import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import pdb

import argparse
parser = argparse.ArgumentParser(prog='./ZLumi')
parser.add_argument(
    '-o', '--output', default='./ZLumi',
    help='specify output dir'
)
args = parser.parse_args()

output = args.output

if os.path.isdir(output):
    print("output dir already exists, please remove or specify another name")
    exit()
os.mkdir(output)

zcountsH = pd.read_csv('/nfs/dust/cms/user/dwalter/CMSSW_10_2_13/src/ZCounting/ZHarvester/ZCounts2017_Pt30_v1/H/result.csv')
zcountsB = pd.read_csv('/nfs/dust/cms/user/dwalter/CMSSW_10_2_13/src/ZCounting/ZHarvester/ZCounts2017_Pt30_v1/B/result.csv')
zcountsC = pd.read_csv('/nfs/dust/cms/user/dwalter/CMSSW_10_2_13/src/ZCounting/ZHarvester/ZCounts2017_Pt30_v1/C/result.csv')
zcountsD = pd.read_csv('/nfs/dust/cms/user/dwalter/CMSSW_10_2_13/src/ZCounting/ZHarvester/ZCounts2017_Pt30_v1/D/result.csv')
zcountsE = pd.read_csv('/nfs/dust/cms/user/dwalter/CMSSW_10_2_13/src/ZCounting/ZHarvester/ZCounts2017_Pt30_v1/E/result.csv')
zcountsF = pd.read_csv('/nfs/dust/cms/user/dwalter/CMSSW_10_2_13/src/ZCounting/ZHarvester/ZCounts2017_Pt30_v1/F/result.csv')

labels = ['(shape, shape)', '(shape, analytic)', '(analytic, shape)', '(analytic, analytic)']
colors = ['']

linewidth = 1.5
rcParams['axes.linewidth'] = linewidth
rcParams['xtick.major.width'] = linewidth
rcParams['ytick.major.width'] = linewidth
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']


#plt.xkcd(0, 0, 0)

### control plots on efficiencies
for eff in ('effGloB', 'effGloE',
            'effSelB', 'effSelE',
            'effHLTB', 'effHLTE',):
    plt.clf()
    fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    fig.subplots_adjust(hspace=0.05, bottom=0.2, top=0.95)

    for zcounts, era, col in ((zcountsB, 'B', 'blue'),
                              (zcountsC, 'C', 'orange'),
                              (zcountsD, 'D', 'green'),
                              (zcountsE, 'E', 'red'),
                              (zcountsF, 'F', 'purple'),
                              (zcountsH, 'H', 'cyan'),
                              ):

        ax[0].plot(zcounts[eff], label='2017{0}'.format(era), marker='o', linewidth=2, color=col)
        ax[1].plot(zcounts[eff]/zcounts[eff][0], marker='o', linewidth=2, color=col)

    ax[0].set_ylabel(eff[:-1])
    yWidth = (ax[0].get_ylim()[1] - ax[0].get_ylim()[0])
    ax[0].set(xlim=(-0.5, 3.5), ylim=(ax[0].get_ylim()[0] - yWidth*0.125, ax[0].get_ylim()[1] + yWidth*0.125))

    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width * 0.85, box.height])

    leg = ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    leg.get_frame().set_linewidth(linewidth)
    leg.get_frame().set_edgecolor("black")

    if eff[-1] == 'E':
        ax[0].text(0.05, 0.1, 'endcap', transform=ax[0].transAxes)
    elif eff[-1] == 'B':
        ax[0].text(0.05, 0.1, 'barrel', transform=ax[0].transAxes)

    box = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0, box.width * 1.1, box.height])

    ax[1].set_xticks(range(4) )
    ax[1].set_xticklabels(labels, rotation=30)
    ax[1].set_ylabel("<"+eff[:-1]+">")
    yWidth = (ax[1].get_ylim()[1] - ax[1].get_ylim()[0])
    ax[1].set(ylim=(ax[1].get_ylim()[0] - yWidth * 0.125, ax[1].get_ylim()[1] + yWidth * 0.125))


    print("create {0}".format(output+"/plots_{0}.png".format(eff)))
    plt.savefig(output+"/plots_{0}.png".format(eff))


### control plots on chi2 values
for eff in ('effGloB_chi2Pass', 'effGloE_chi2Pass',
            'effSelB_chi2Pass', 'effSelE_chi2Pass',
            'effHLTB_chi2Pass', 'effHLTE_chi2Pass',
            'effGloB_chi2Fail', 'effGloE_chi2Fail',
            'effSelB_chi2Fail', 'effSelE_chi2Fail',
            'effHLTB_chi2Fail', 'effHLTE_chi2Fail',
            ):
    plt.clf()
    fig, ax = plt.subplots(1, 1)
    fig.subplots_adjust(bottom=0.2,)

    for zcounts, era, col in ((zcountsB, 'B', 'blue'),
                              (zcountsC, 'C', 'orange'),
                              (zcountsD, 'D', 'green'),
                              (zcountsE, 'E', 'red'),
                              (zcountsF, 'F', 'purple'),
                              (zcountsH, 'H', 'cyan'),
                              ):

        ax.plot(zcounts[eff], label='2017{0}'.format(era), marker='o', linewidth=2, color=col)

    ax.set_ylabel("chi2")
    yWidth = (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.set(xlim=(-0.5, 3.5), ylim=(ax.get_ylim()[0] - yWidth*0.125, ax.get_ylim()[1] + yWidth*0.125))

    leg = ax.legend(loc='center', bbox_to_anchor=(0.95, 0.8))

    leg.get_frame().set_linewidth(linewidth)
    leg.get_frame().set_edgecolor("black")

    ax.text(0.05, 0.9, '{0} - {1} - {2}'.format(eff[3:6], 'barrel' if eff[6] == 'B' else 'endcap', eff[-4:]), transform=ax.transAxes)

    ax.set_xticks(range(4) )
    ax.set_xticklabels(labels, rotation=30)

    print("create {0}".format(output+"/plots_{0}.png".format(eff)))
    plt.savefig(output+"/plots_{0}.png".format(eff))


#brilcalc lumi -c web -i /eos/home-d/dwalter/CMSSW_10_6_4/src/ZCounting/TnPPairTreeProducer/production/Cert_306896-307082_13TeV_PromptReco_Collisions17_JSON_LowPU.txt --hltpath="HLT_HIMu17_v*" --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb
LumiH = 199.8
LumiH_Erel = 0.017

Lumi = {'B': 4794, 'C': 9632, 'D': 4248, 'E': 9315, 'F': 13540}
LumiRun2 = 41528

ZLumiRun2 = 0.
ZLumiRun2_Etot = 0.
for zcounts, era in ((zcountsB, 'B'), (zcountsC, 'C'), (zcountsD, 'D'), (zcountsE, 'E'), (zcountsF, 'F')):





    ### luminosity
    ZLumi = LumiH * zcounts['NZDelivMC'] / zcountsH['NZDelivMC']

    print("Lumi({}) shape variations: ".format(era))
    print(ZLumi)

    ZLumi_EstatH = ZLumi[0] * float(zcountsH['NZDelivMC_EStat'][0])/zcountsH['NZDelivMC'][0]
    ZLumi_EStat = ZLumi[0] * float(zcounts['NZDelivMC_EStat'][0])/zcounts['NZDelivMC'][0]
    ZLumi_ELumH = ZLumi[0] * LumiH_Erel
    ZLumi_EShapes = 0.005 * ZLumi[0]
    ZLumi_Etot = np.sqrt(ZLumi_ELumH**2 + ZLumi_EShapes**2 + ZLumi_EstatH**2 + ZLumi_EStat)
    print("ZLumi({0}) = {1} +- {2} (lumi) +- {3} (shapes) {4} (statLowPU) +- {5} (statHighPU)".format(era, ZLumi[0], ZLumi_ELumH, ZLumi_EShapes, ZLumi_EstatH, ZLumi_EStat))
    print("ZLumi({0}) = {1} +- {2}".format(era, ZLumi[0], ZLumi_Etot))
    print("Lumi({0}) = {1} +- {2}".format(era, Lumi[era], Lumi[era]*0.023))

    ZLumiRun2 += ZLumi[0]
    ZLumiRun2_Etot += ZLumi_Etot

print("ZLumi(run2) = {0} +- {1}".format(ZLumiRun2, ZLumiRun2_Etot))
print("Lumi(run2) = {0} +- {1}".format(LumiRun2, LumiRun2*0.023))