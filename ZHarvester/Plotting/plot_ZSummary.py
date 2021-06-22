import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
import uncertainties as unc
import pdb
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime, unorm

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("--rates", required=True, nargs='+', help="Nominator csv file with z rates per measurement")
parser.add_argument("--xsec",  required=True, type=str,
    help="csv file with z rates per measurement where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

# --- settings
secPerLS=float(23.3)
labelsize = 12.5
textsize = 15

# --- uncertainties on PHYSICS luminosity
include_unc_PHYSICS = True
unc_2016 = np.sqrt((0.012)**2 + (0.017)**2 - 2*0.012*0.017*0.26)
unc_2017 = np.sqrt((0.023)**2 + (0.017)**2 - 2*0.023*0.017*0.76)
unc_2018 = np.sqrt((0.025)**2 + (0.017)**2 - 2*0.025*0.017*0.43)

print("relative uncertainties attributed to PHYSICS: ")
print("2016: "+str(unc_2016))
print("2017: "+str(unc_2017))
print("2018: "+str(unc_2018))

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Lato', 'DejaVu Sans',
                               'Lucida Grande', 'Verdana']

########## Data Acquisition ##########

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

data['ZBBeff_mc'] = data['ZBBeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['ZBEeff_mc'] = data['ZBEeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['ZEEeff_mc'] = data['ZEEeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))

data['zDelBB_mc'] = data['zDelBB_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDelBE_mc'] = data['zDelBE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDelEE_mc'] = data['zDelEE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDel_mc'] = data['zDelBB_mc'] + data['zDelBE_mc'] + data['zDelEE_mc']

data['zLumiBB_mc'] = data['zDelBB_mc'] / xsecBB
data['zLumiBE_mc'] = data['zDelBE_mc'] / xsecBE
data['zLumiEE_mc'] = data['zDelEE_mc'] / xsecEE
data['zLumi_mc'] = data['zDel_mc'] / (xsecBB+xsecBE+xsecEE)

data['zDel_mc'] = data['zDel_mc'].apply(lambda x: x.n)
data['zLumi_mc'] = data['zLumi_mc'].apply(lambda x: x.n)
data['zLumiBB_mc'] = data['zLumiBB_mc'].apply(lambda x: x.n)
data['zLumiBE_mc'] = data['zLumiBE_mc'].apply(lambda x: x.n)
data['zLumiEE_mc'] = data['zLumiEE_mc'].apply(lambda x: x.n)

data['time'] = (data['tdate_begin']+data['tdate_end'])//2

# data = data[['run', 'fill','zLumi_mc','zLumiBB_mc','zLumiBE_mc','zLumiEE_mc','time','lumiRec']]

data = data[data['lumiRec'] > 0.]
data = data[data['zLumi_mc'] > 0.]

data['zLumi_mc_to_dLRec'] = data['zLumi_mc'] / data['lumiRec']
data['zLumiBB_mc_to_dLRec'] = data['zLumiBB_mc'] / data['lumiRec']
data['zLumiBE_mc_to_dLRec'] = data['zLumiBE_mc'] / data['lumiRec']
data['zLumiEE_mc_to_dLRec'] = data['zLumiEE_mc'] / data['lumiRec']
data['weightLumi'] = data['lumiRec']

print("analyze {0} fb^-1 of data".format(data['weightLumi'].sum()/1000.))

def make_hist(
    df,
    run_range=None,
    lumi_name='zLumi_mc_to_dLRec',
    sumN=1,
    label="ZCount / PHYSICS",
    saveas="zcount",
    title=None
):

    # if sumN == 1:
    #     label += " by measurement".format(sumN)
    # else:
    #     label += " by {0} measurements".format(sumN)


    lefttitle = "$\sqrt{s}=13\,\mathrm{TeV}$"

    if title:
        lefttitle += " $(\mathrm{"+title+"})$"

    saveas = str(sumN) + "_" + saveas

    if run_range:
        data = df.loc[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])]
        if len(df) ==0:
            return
    else:
        data = df

    data = data.sort_values(['run','time'])

    # --- sum up each sumN rows
    lumiratio = data[lumi_name].groupby(data.index // sumN).sum()/sumN
    ZBBeff_mc = data["ZBBeff_mc"].groupby(data.index // sumN).sum()/sumN
    ZBEeff_mc = data["ZBEeff_mc"].groupby(data.index // sumN).sum()/sumN
    ZEEeff_mc = data["ZEEeff_mc"].groupby(data.index // sumN).sum()/sumN
    time = data['time'].groupby(data.index // sumN).mean()
    run = data['run'].groupby(data.index // sumN).mean()
    fill = data['fill'].groupby(data.index // sumN).mean()
    weight = data['weightLumi'].groupby(data.index // sumN).sum()

    # --- make histogram
    # mean and std without outliers
    mean = lumiratio[(lumiratio<2.0) & (lumiratio>0.5)].mean()
    std = lumiratio[(lumiratio<2.0) & (lumiratio>0.5)].std()
    range = (mean-2*std,mean+2*std)

    for weighted in (True, False):
        plt.clf()
        fig = plt.figure()
        fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.125)
        ax = fig.add_subplot(111)
        if weighted:
            nEntries, bins, _ = ax.hist(lumiratio.values, weights=weight.values, bins=50, range=(mean-2*std,mean+2*std))
            ax.set_ylabel("Integrated luminosity [pb$^{-1}$]", fontsize=textsize)
        else:
            nEntries, bins, _ = ax.hist(lumiratio.values, bins=50, range=(mean-2*std,mean+2*std))
            ax.set_ylabel("Number of entries", fontsize=textsize)

        ax.set_xlabel(label, fontsize=textsize)
        ax.text(0.04, 0.97, "CMS", verticalalignment='top', transform=ax.transAxes, weight="bold", fontsize=textsize)
        ax.text(0.15, 0.97, "Work in progress", verticalalignment='top', transform=ax.transAxes,style='italic', fontsize=textsize)
        ax.text(0.04, 0.91, lefttitle, verticalalignment='top', transform=ax.transAxes,style='italic', fontsize=textsize)
        ax.text(0.7, 0.97, "mean = {0}".format(round(mean,3)), verticalalignment='top', transform=ax.transAxes, fontsize=textsize)
        ax.text(0.7, 0.91, "std = {0}".format(round(std,3)), verticalalignment='top', transform=ax.transAxes, fontsize=textsize)

        histname = "/hist_weighted_"+saveas if weighted else "/hist_"+saveas
        print("save histogram as {0}".format(outDir+histname))
        plt.xticks(fontsize = labelsize)
        plt.yticks(fontsize = labelsize)

        plt.savefig(outDir+histname)
        plt.close()

    # --- make scatter
    rangey = (mean-4*std,mean+4*std)
    for xx, xlabel, suffix1 in (
        # (run.values, "Run number", "run"),
        # (fill.values, "Fill number", "fill"),
        # # (time.values, "Time", "time"),
        (weight.cumsum().values/1000., "Integrated PHYSICS luminosity [fb$^{-1}$]", "lumi"),
    ):
        rangex = min(xx)-(max(xx)-min(xx))*0.01, max(xx)+(max(xx)-min(xx))*0.01

        for yy, ylabel, suffix in (
            (lumiratio.values, label, "lumi"),
            # (np.array([x.n for x in ZBBeff_mc.values]), "Z efficiency (BB)", "ZeffBB"),
            # (np.array([x.n for x in ZBEeff_mc.values]), "Z efficiency (BE)", "ZeffBE"),
            # (np.array([x.n for x in ZEEeff_mc.values]), "Z efficiency (EE)", "ZeffEE"),
        ):
            mean = np.mean(yy)
            std = np.std(yy)
            rangey = (mean-4*std,mean+4*std)

            plt.clf()
            fig = plt.figure(figsize=(10.0,4.0))
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.08, right=0.99, top=0.99, bottom=0.15)

            # plot uncertinty bar attributed to PHYSICS luminosity
            if suffix1 in ("lumi", ) and include_unc_PHYSICS:

                starts = np.array([rangex[0],])
                widths = np.array([abs(rangex[1] - rangex[0]),])

                if title and "2016" in title:
                    heights = np.array([unc_2016*2,])
                    bottoms = np.array([1. - unc_2016,])
                elif title and "2017" in title:
                    heights = np.array([unc_2017*2,])
                    bottoms = np.array([1. - unc_2017,])
                elif title and "2018" in title:
                    heights = np.array([unc_2018*2,])
                    bottoms = np.array([1. - unc_2018,])
                else:
                    starts = np.array([rangex[0], 36.33, 78.19])
                    heights = np.array([unc_2016*2, unc_2017*2, unc_2018*2])
                    widths = np.array([abs(rangex[0])+36.33, 41.86, rangex[1]-78.19])
                    bottoms = np.array([1. - unc_2016, 1. - unc_2017, 1. - unc_2018])

                ax.bar(starts, height=heights, width=widths, bottom=bottoms, align='edge',
                    color="grey", alpha=0.4, hatch='/', zorder=4, #, alpha=0.6
                    label="PHYSICS uncertainty")

            ss = 20 * weight.values / np.mean(weight.values)
            ax.scatter(xx, yy, s=ss, marker='.', color='green', zorder=2, label="Measurement")

            ax.text(0.03, 0.97, "CMS", verticalalignment='top', transform=ax.transAxes, weight="bold", fontsize=textsize)
            ax.text(0.10, 0.97, "Work in progress", verticalalignment='top', transform=ax.transAxes,style='italic', fontsize=textsize)
            ax.text(0.98, 0.97, lefttitle, verticalalignment='top', transform=ax.transAxes,style='italic', fontsize=textsize, horizontalalignment='right')

            ax.set_xlabel(xlabel, fontsize=textsize)
            ax.set_ylabel(ylabel, fontsize=textsize)
            ax.set_ylim(rangey)
            ax.set_xlim(rangex)
            if suffix1 in ("lumi", ):
                # plot horizontal line at 1
                ax.plot(np.array(rangex), np.array([1.,1.]), 'k--', linewidth=1, zorder=3)

            print("save scatter as {0}".format(outDir+"/scatter_"+suffix+"_"+saveas))
            plt.xticks(fontsize = labelsize)
            plt.yticks(fontsize = labelsize)
            plt.legend(loc='lower right', ncol=2, markerscale=3, scatteryoffsets=[0.5], fontsize=textsize, frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
            plt.savefig(outDir+"/scatter_"+suffix+"_"+saveas)
            plt.close()


make_hist(data, lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", saveas="zcount.pdf", title="Run\ II")
make_hist(data, lumi_name='zLumiBB_mc_to_dLRec', label="ZCount(BB) / PHYSICS", saveas="zcountBB.pdf")
make_hist(data, lumi_name='zLumiBE_mc_to_dLRec', label="ZCount(BE) / PHYSICS", saveas="zcountBE.pdf")
make_hist(data, lumi_name='zLumiEE_mc_to_dLRec', label="ZCount(EE) / PHYSICS", saveas="zcountEE.pdf")

make_hist(data, run_range=(272007,294645), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", saveas="2016_zcount.pdf", title="2016")
make_hist(data, run_range=(272007,278769), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", saveas="2016preVFP_zcount.pdf", title="2016\ pre\ VFP")
make_hist(data, run_range=(278769,294645), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", saveas="2016postVFP_zcount.pdf", title="2016\ post\ VFP")
make_hist(data, run_range=(297046,306462), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", saveas="2017_zcount.pdf", title="2017")
make_hist(data, run_range=(315252,325175), lumi_name='zLumi_mc_to_dLRec', label="ZCount / PHYSICS", saveas="2018_zcount.pdf", title="2018")
