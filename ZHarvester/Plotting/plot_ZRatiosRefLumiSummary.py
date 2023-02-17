import os,sys
import ROOT
import argparse
import pandas as pd
import numpy as np
import uncertainties as unc
import pdb
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import unorm

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

from common import parsing
from common.logging import child_logger
log = child_logger(__name__)

parser = parsing.parser_plot()

parser.add_argument("--rates", required=True, nargs='+', help="Nominator csv file with z rates perLS")
parser.add_argument("--xsec",  required=True, type=str,
    help="csv file with z rates per measurement where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("--refLumi",  required=True, nargs='+', help="give a ByLs.csv as input for reference Luminosity")
args = parser.parse_args()

outDir = args.output
if not os.path.isdir(outDir):
    os.mkdir(outDir)

# --- settings

secPerLS=float(23.3)
textsize = 20


########## Data Acquisition ##########

# --- reference luminosity
print("convert reference data from ByLS.csv")
data_ref = []
for csv in args.refLumi:
    lumiFile=open(csv)
    lumiLines=lumiFile.readlines()
    data_ref.append(pd.read_csv(csv, sep=',',low_memory=False,
        skiprows=lambda x: lumiLines[x].startswith('#') and not lumiLines[x].startswith('#run')
        ))
data_ref = pd.concat(data_ref, ignore_index=True, sort=False)
if 'recorded(/ub)' in data_ref.columns.tolist():      #convert to /pb
    data_ref['recorded(/ub)'] = data_ref['recorded(/ub)'].apply(lambda x:x / 1000000.)
    data_ref = data_ref.rename(index=str, columns={'recorded(/ub)':'recorded(/pb)' })
elif 'recorded(/fb)' in data_ref.columns.tolist():      #convert to /pb
    data_ref['recorded(/fb)'] = data_ref['recorded(/fb)'].apply(lambda x:x * 1000.)
    data_ref = data_ref.rename(index=str, columns={'recorded(/fb)':'recorded(/pb)' })

data_ref['fill'] = pd.to_numeric(data_ref['#run:fill'].str.split(':',expand=True)[1])
data_ref['run'] = pd.to_numeric(data_ref['#run:fill'].str.split(':',expand=True)[0])
data_ref['ls'] = pd.to_numeric(data_ref['ls'].str.split(':',expand=True)[0])
data_ref = data_ref.drop(['#run:fill','hltpath','source'],axis=1)
data_ref = data_ref.sort_values(['fill','run','ls','recorded(/pb)'])
data_ref = data_ref.drop_duplicates(['fill','run','ls'])

# currentYear=2017
# data_ref['time'] = data_ref['time'].apply(lambda x: to_RootTime(x,currentYear)).astype(float)
# data_ref['timeE'] = data_ref['time'] - data_ref['time']

data_ref['dLRec(/pb)'] = data_ref['recorded(/pb)']
data_ref = data_ref[['run','ls','dLRec(/pb)']]		#Keep only what you need

# --- get Z xsec
print("get Z cross section")
data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])

data_xsec['yieldBB'] = data_xsec['yieldBB'].apply(lambda x: unorm(x))
data_xsec['yieldBE'] = data_xsec['yieldBE'].apply(lambda x: unorm(x))
data_xsec['yieldEE'] = data_xsec['yieldEE'].apply(lambda x: unorm(x))

data_xsec['ZBBeff_mc'] = data_xsec['ZBBeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data_xsec['ZBEeff_mc'] = data_xsec['ZBEeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data_xsec['ZEEeff_mc'] = data_xsec['ZEEeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))

data_xsec['zYieldBB_purity'] = data_xsec['zYieldBB_purity'].apply(lambda x: unc.ufloat_fromstr(x))
data_xsec['zYieldBE_purity'] = data_xsec['zYieldBE_purity'].apply(lambda x: unc.ufloat_fromstr(x))
data_xsec['zYieldEE_purity'] = data_xsec['zYieldEE_purity'].apply(lambda x: unc.ufloat_fromstr(x))

data_xsec['zYieldBB_purity'] = data_xsec['zYieldBB_purity'].apply(lambda x: unc.ufloat(1,0.) if x.n == 1.0 else x)
data_xsec['zYieldBE_purity'] = data_xsec['zYieldBE_purity'].apply(lambda x: unc.ufloat(1,0.) if x.n == 1.0 else x)
data_xsec['zYieldEE_purity'] = data_xsec['zYieldEE_purity'].apply(lambda x: unc.ufloat(1,0.) if x.n == 1.0 else x)

# delivered Z rate with applying purity
data_xsec['zDelBB_mc'] = data_xsec['yieldBB'] / data_xsec['ZBBeff_mc'] * data_xsec['zYieldBB_purity']
data_xsec['zDelBE_mc'] = data_xsec['yieldBE'] / data_xsec['ZBEeff_mc'] * data_xsec['zYieldBE_purity']
data_xsec['zDelEE_mc'] = data_xsec['yieldEE'] / data_xsec['ZEEeff_mc'] * data_xsec['zYieldEE_purity']
data_xsec['zDel_mc'] = data_xsec['zDelBB_mc'] + data_xsec['zDelBE_mc'] + data_xsec['zDelEE_mc']

xsecBB = sum(data_xsec['zDelBB_mc'])/sum(data_xsec['lumiRec'])
xsecBE = sum(data_xsec['zDelBE_mc'])/sum(data_xsec['lumiRec'])
xsecEE = sum(data_xsec['zDelEE_mc'])/sum(data_xsec['lumiRec'])

# --- z luminosity
print("get Z luminosity")
data = pd.concat([pd.read_csv(csv, sep=',',low_memory=False) for csv in args.rates], ignore_index=True, sort=False)

# data['yieldBB'] = data['yieldBB'].apply(lambda x: unorm(x))
# data['yieldBE'] = data['yieldBE'].apply(lambda x: unorm(x))
# data['yieldEE'] = data['yieldEE'].apply(lambda x: unorm(x))
#
# data['effBB_mc'] = data['effBB_mc'].apply(lambda x: unc.ufloat_fromstr(x))
# data['effBE_mc'] = data['effBE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
# data['effEE_mc'] = data['effEE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
#
# data['zYieldBB_purity'] = data['zYieldBB_purity'].apply(lambda x: unc.ufloat_fromstr(x))
# data['zYieldBE_purity'] = data['zYieldBE_purity'].apply(lambda x: unc.ufloat_fromstr(x))
# data['zYieldEE_purity'] = data['zYieldEE_purity'].apply(lambda x: unc.ufloat_fromstr(x))
#
# data['zYieldBB_purity'] = data['zYieldBB_purity'].apply(lambda x: unc.ufloat(1,0.) if x.n == 1.0 else x)
# data['zYieldBE_purity'] = data['zYieldBE_purity'].apply(lambda x: unc.ufloat(1,0.) if x.n == 1.0 else x)
# data['zYieldEE_purity'] = data['zYieldEE_purity'].apply(lambda x: unc.ufloat(1,0.) if x.n == 1.0 else x)
#
# # delivered Z rate with applying purity
# data['zDelBB_mc'] = data['yieldBB'] / data['effBB_mc'] * data['zYieldBB_purity']
# data['zDelBE_mc'] = data['yieldBE'] / data['effBE_mc'] * data['zYieldBE_purity']
# data['zDelEE_mc'] = data['yieldEE'] / data['effEE_mc'] * data['zYieldEE_purity']
# data['zDel_mc'] = data['zDelBB_mc'] + data['zDelBE_mc'] + data['zDelEE_mc']

data['zDelBB_mc'] = data['zDelBB_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDelBE_mc'] = data['zDelBE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDelEE_mc'] = data['zDelEE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDel_mc'] = data['zDelBB_mc'] + data['zDelBE_mc'] + data['zDelEE_mc']

data['zLumiBB_mc'] = data['zDelBB_mc'] / xsecBB
data['zLumiBE_mc'] = data['zDelBE_mc'] / xsecBE
data['zLumiEE_mc'] = data['zDelEE_mc'] / xsecEE
data['zLumi_mc'] = data['zDel_mc'] / (xsecBB+xsecBE+xsecEE)

# sort out lumi section withou any counts
data['zDel_mc'] = data['zDel_mc'].apply(lambda x: x.n)
data['zLumi_mc'] = data['zLumi_mc'].apply(lambda x: x.n)
data['zLumiBB_mc'] = data['zLumiBB_mc'].apply(lambda x: x.n)
data['zLumiBE_mc'] = data['zLumiBE_mc'].apply(lambda x: x.n)
data['zLumiEE_mc'] = data['zLumiEE_mc'].apply(lambda x: x.n)

data = data[['run','ls','zLumi_mc','zLumiBB_mc','zLumiBE_mc','zLumiEE_mc','time']]

# merge data_ref and data; keep run,ls that are in both
print("merge Z luminosity and reference luminosity")
merged = pd.merge(data_ref, data, on=['run','ls'])

merged = merged[merged['dLRec(/pb)'] > 0]
merged = merged[merged['zLumi_mc'] > 0]

merged['zLumi_mc_to_dLRec'] = merged['zLumi_mc'] / merged['dLRec(/pb)']
merged['zLumiBB_mc_to_dLRec'] = merged['zLumiBB_mc'] / merged['dLRec(/pb)']
merged['zLumiBE_mc_to_dLRec'] = merged['zLumiBE_mc'] / merged['dLRec(/pb)']
merged['zLumiEE_mc_to_dLRec'] = merged['zLumiEE_mc'] / merged['dLRec(/pb)']
# merged['weightLumi'] = (merged['zLumiInst_mc']+merged['dLRec(/pb)'])/2.
merged['weightLumi'] = merged['dLRec(/pb)']

print("analyze {0} fb^-1 of data".format(merged['weightLumi'].sum()/1000.))

def make_hist(
    df,
    run_range=None,
    lumi_name='zLumi_mc_to_dLRec',
    sumN=200,
    label="ZCount / PHYSICS",
    saveas="zcount.png"
):

    label += " by {0} ls".format(sumN)

    saveas = str(sumN) + "_" + saveas

    if run_range:
        data = df.loc[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])]
        if len(df) ==0:
            return
    else:
        data = df

    data = data.sort_values(['run','ls'])

    # --- sum up each sumN rows
    res = data[lumi_name].groupby(data.index // sumN).sum()/sumN
    time = data['time'].groupby(data.index // sumN).mean()
    run = data['run'].groupby(data.index // sumN).mean()
    weight = data['weightLumi'].groupby(data.index // sumN).sum()

    # --- make histogram
    # mean and std without outliers
    mean = res[(res<2.0) & (res>0.5)].mean()
    std = res[(res<2.0) & (res>0.5)].std()
    range = (mean-2*std,mean+2*std)

    for weighted in (True, False):
        plt.clf()
        fig = plt.figure()
        fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.1)
        ax = fig.add_subplot(111)
        if weighted:
            nEntries, bins, _ = ax.hist(res.values, weights=weight.values, bins=50, range=(mean-2*std,mean+2*std))
            ax.set_ylabel("Integrated luminosity [pb$^{-1}$]", fontsize=textsize)
        else:
            nEntries, bins, _ = ax.hist(res.values, bins=50, range=(mean-2*std,mean+2*std))
            ax.set_ylabel("Number of entries", fontsize=textsize)

        ax.set_xlabel(label, fontsize=textsize)
        ax.text(0.04, 0.97, "CMS", verticalalignment='top', transform=ax.transAxes, weight="bold", fontsize=textsize)
        ax.text(0.15, 0.97, "Work in progress", verticalalignment='top', transform=ax.transAxes,style='italic', fontsize=textsize)
        ax.text(0.8, 0.97, "mean = {0}".format(round(mean,3)), verticalalignment='top', transform=ax.transAxes, fontsize=textsize)
        ax.text(0.8, 0.92, "std = {0}".format(round(std,3)), verticalalignment='top', transform=ax.transAxes, fontsize=textsize)

        histname = "/hist_weighted_"+saveas if weighted else "/hist_"+saveas
        print("save histogram as {0}".format(outDir+histname))
        plt.savefig(outDir+histname)
        plt.close()

    # --- make scatter
    rangey = (mean-4*std,mean+4*std)
    for xx, xlabel, suffix in (
        (run.values,"Run number", "run"),
        # (time.values, "Time", "time"),
        (weight.cumsum().values/1000., "Integrated PHYSICS luminosity [fb$^{-1}$]", "lumi")
    ):
        rangex = min(xx)-(max(xx)-min(xx))*0.01, max(xx)+(max(xx)-min(xx))*0.01

        plt.clf()
        fig = plt.figure(figsize=(10.0,4.0))
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.08, right=0.99, top=0.99, bottom=0.15)
        ax.scatter(xx, res.values, marker='.', color='green')

        ax.text(0.03, 0.97, "CMS", verticalalignment='top', transform=ax.transAxes, weight="bold", fontsize=textsize)
        ax.text(0.10, 0.97, "Work in progress", verticalalignment='top', transform=ax.transAxes,style='italic', fontsize=textsize)

        ax.set_xlabel(xlabel, fontsize=textsize)
        ax.set_ylabel(label, fontsize=textsize)
        ax.set_ylim(rangey)
        ax.set_xlim(rangex)
        ax.plot(np.array(rangex), np.array([1.,1.]), 'k--', linewidth=1)
        print("save scatter as {0}".format(outDir+"/scatter_"+suffix+"_"+saveas))
        plt.savefig(outDir+"/scatter_"+suffix+"_"+saveas)
        plt.close()


make_hist(merged, lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS",sumN=100, saveas="zcount.png")
make_hist(merged, lumi_name='zLumiBB_mc_to_dLRec', label="ZCount(BB) / PHYSICS",sumN=100, saveas="zcountBB.png")
make_hist(merged, lumi_name='zLumiBE_mc_to_dLRec', label="ZCount(BE) / PHYSICS",sumN=100, saveas="zcountBE.png")
make_hist(merged, lumi_name='zLumiEE_mc_to_dLRec', label="ZCount(EE) / PHYSICS",sumN=100, saveas="zcountEE.png")

make_hist(merged, run_range=(272007,294645), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", sumN=100, saveas="2016_zcount.png")
make_hist(merged, run_range=(272007,278769), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", sumN=100, saveas="2016preVFP_zcount.png")
make_hist(merged, run_range=(278769,294645), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", sumN=100, saveas="2016postVFP_zcount.png")
make_hist(merged, run_range=(297046,306462), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", sumN=100, saveas="2017_zcount.png")
make_hist(merged, run_range=(315252,325175), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", sumN=100, saveas="2018_zcount.png")
