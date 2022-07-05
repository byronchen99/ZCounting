import ROOT
import pandas as pd
import numpy as np
import argparse
import pdb
import os, sys
import math
import uncertainties as unc
from scipy.optimize import curve_fit
from scipy import stats

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

sys.path.append(os.getcwd())
print(os.getcwd())

from python.corrections import apply_muon_prefire, apply_ECAL_prefire, apply_pileup_correction

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime, cms, preliminary, text, workinprogress, unorm, linear

pd.options.mode.chained_assignment = None

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)

textsize = 16
markersize = 4.0

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino",],
    "font.size": textsize,
    'text.latex.preamble': [r"""\usepackage{bm}"""]
})

mpl.rcParams.update({
    "legend.fontsize" : "medium",
    "axes.labelsize" : "medium",
    "axes.titlesize" : "medium",
    "xtick.labelsize" : "medium",
    "ytick.labelsize" : "medium",
})
colors = ["#e74c3c","#f1c40f","#16a085","#27ae60","#2980b9","#8e44ad"]

parser = argparse.ArgumentParser()

parser.add_argument("-r","--rates", required=True, type=str, help="csv file with z rates per Measurement")
parser.add_argument("-x","--xsec",  type=str, help="csv file where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")

args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

# ------------------------------------------------------------------------------
def make_plots(df,
    yAxis,
    yLabel="sigma",
    xAxis='lumi',
    run_range=None,
    title="",
    year="2017",
    region='inclusive',
    plot_fills=False,
    plot_all=False,
    nround=0,
    normalized=False,
    resource=""
):
    """
    valid xAxis: 'lumi', 'pileUp', 'measurement', 'time'
    """

    if run_range:
        data = df.loc[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])]
        if len(df) ==0:
            return
    else:
        data = df

    data['y_Err'] = data[yAxis].apply(lambda x: x.s)
    data['y'] = data[yAxis].apply(lambda x: x.n)
    if yAxis.replace("_mc","") in data.keys():
        data['y0_Err'] = data[yAxis.replace("_mc","")].apply(lambda x: x.s)
        data['y0'] = data[yAxis.replace("_mc","")].apply(lambda x: x.n)
    else:
        data['y0_Err'] = data['y_Err']
        data['y0'] = data['y']

    if sum(data['y'].isnull()) > 0:
        print(">>> sort out {0} points with nan".format(sum(data['y'].isnull())))
        data = data.loc[~data['y'].isnull()]

    data['relstat'] = data['y_Err'] / data['y']

    # x_step: intervall in which the measurements are collected in one point

    if xAxis == 'lumi':
        xAxis = 'lumiRecInst'
        data[xAxis] = data['lumiRec'] / data['timewindow']
        xTitle="Inst. luminosity $[\\mathrm{nb}^{-1}\\mathrm{s}^{-1}]$"
        x_step = 0.2
    elif xAxis == 'time':
        data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2
        #data['time'] -= 7200
        data['time'] -= (ROOT.TDatime(2042,01,01,12,00,00).Convert() - ROOT.TDatime(2017,01,01,12,00,00).Convert())
        xTitle="time"
        x_step = 1
    elif xAxis == 'pileUp':
        xTitle="average number of pileup interactions"
        x_step = 2
    else:
        xTitle=xAxis
        x_step = 0.2

    if yLabel == "sigma":
        yLabel = "$\\sigma^{\\mathrm{fid}}_{\\mathrm{Z}}$\,[pb]"

        if normalized:
            data[yAxis] /= np.mean(data[yAxis].values)
            data[yAxis.replace("_mc","")] /= np.mean(data[yAxis.replace("_mc","")].values)
            yLabel = "$\\sigma^{\\mathrm{fid}}_{\\mathrm{Z}} / \\langle \\sigma^{\\mathrm{fid}}_{\\mathrm{Z}} \\rangle$ "

    elif resource != "":
        # plot an efficiency on y axis - load necessery resources

        # >>>  load MC info and format
        if len(year.split(" ")) == 1:
            suffix = year
        elif len(year.split(" ")) == 0:
            print("Invalid year")
        elif year.split(" ")[0] == "2016":
            if year == "2016 pre VFP" or year.split(" ")[1] in ("B","C","D","E","F"):
                suffix="2016preVFP"
            elif year == "2016 post VFP" or year.split(" ")[1] in ("G","H"):
                suffix="2016postVFP"
        else:
            suffix = year.split(" ")[0]

        basepath = "/nfs/dust/cms/user/dwalter/CMSSW_10_2_20_UL/src/potato-zcount/plots/"
        genInfoFile = basepath+'GenInfo-V13_02-d20211012-t180358/infoMC_gen_{0}.json'.format(suffix)
        recoInfoFile = basepath+'RecoInfo-V13_02-d20211012-t180324/infoMC_reco_{0}.json'.format(suffix)

        infoReco = pd.read_json(recoInfoFile,orient="index")
        infoTrue = pd.read_json(genInfoFile,orient="index")
        infoReco['bin']= infoReco.index
        infoTrue['bin']= infoTrue.index

        infoReco["True"] = False
        infoTrue["True"] = True

        info = pd.merge(infoTrue, infoReco, how='outer')
        info = info.loc[info['bin'] != 'all']
        info[['lo','hi']] = info['bin'].str.split("To", expand=True)
        info['lo'] = info['lo'].apply(lambda x: int(x.replace("nPU","")))
        info['hi'] = info['hi'].apply(lambda x: int(x.replace("Inf","75")))
        info['center'] = info['lo'] + (info['hi'] - info['lo'])/2.

        info.sort_values('center', inplace=True)

        iReco = info.loc[info['True'] == False]
        iTrue = info.loc[info['True'] == True]

        iReco.reset_index(drop=True)
        iTrue.reset_index(drop=True)

    data['y_Err'] = data[yAxis].apply(lambda x: x.s)
    data['y'] = data[yAxis].apply(lambda x: x.n)
    if yAxis.replace("_mc","") in data.keys():
        print(">>> found uncorrected values")
        data['y0_Err'] = data[yAxis.replace("_mc","")].apply(lambda x: x.s)
        data['y0'] = data[yAxis.replace("_mc","")].apply(lambda x: x.n)
    else:
        data['y0_Err'] = data['y_Err']
        data['y0'] = data['y']

    data['relstat'] = data['y_Err'] / data['y']

    # print(">>> sort out {0} points with low statistic".format(len(data['relstat'] > 0.05)))
    # data = data.loc[data["relstat"] < 0.05]

    if plot_fills:
        print(">>> make plots for each fill")
        # For each Fill
        for fill in data.drop_duplicates('fill')['fill'].values:
            dFill = data.loc[data['fill'] == fill]

            subDir = outDir+"/PlotsXSec_"+str(fill)
            if not os.path.isdir(subDir):
                os.mkdir(subDir)
                
            plt.clf()
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.125, hspace=0.0)

            ax1.set_xlabel(xTitle)
            ax1.set_ylabel(yLabel)
            
            ax1.errorbar(dFill[xAxis].values, dFill['y'].values, yerr=data['y_Err'].values, label="Measurement",
                marker="o", linewidth=0, color=colors[0], ecolor=colors[0], elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
                zorder=1)

            leg = ax1.legend(loc="upper left", ncol=3,
                frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
            leg.get_frame().set_linewidth(0.8)
            
            # ax1.set_ylim([minY-yRange*0.05, maxY + yRange*0.5])
            # ax1.set_xlim([xMin, xMax])

            plt.savefig(subDir+"/"+yAxis+"_vs_"+xAxis+".png")
            plt.close()


    if plot_all:
        print(">>> make plot all measurements")
        if xAxis == 'measurement':
            xPoints = np.arange(0,len(data),1.)
        else:
            xPoints = data[xAxis].values

        plt.clf()
        fig = plt.figure(figsize=(10.0,4.0))
        ax1 = fig.add_subplot(111)
        fig.subplots_adjust(left=0.1, right=0.99, top=0.99, bottom=0.15)

        ax1.set_xlabel(xTitle)
        ax1.set_ylabel(yLabel)

        ax1.errorbar(xPoints, data['y'].values, yerr=data['y_Err'].values, label="Measurements",
            marker="o", linewidth=0, color=colors[0], ecolor=colors[0], elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)

        leg = ax1.legend(loc="upper left", ncol=3,
            frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
        leg.get_frame().set_linewidth(0.8)

        # ax1.set_ylim([minY-yRange*0.05, maxY + yRange*0.5])
        # ax1.set_xlim([xMin, xMax])

        outstring = "{0}/{1}_vs_{2}".format(outDir,yAxis,xAxis)
        if run_range:
            outstring += "_run{0}to{1}".format(*run_range)

        plt.savefig(outstring+".png")
        plt.close()


        print(">>> Producing all cross sections")
        print(">>> the simple average cross section is "+str(sum(data['y'].values)/len(data)))

    # make fit taking into account all points
    #print(">>> make fit")
    #func = linear #pol2 #quad
    #popt, pcov = curve_fit(func, data[xAxis], data['y'], sigma=data['y_Err'], absolute_sigma=True)

    #perr = np.sqrt(np.diag(pcov))
    #params = unc.correlated_values(popt, pcov)

    #f = lambda x: func(x, *params)
    #xMC = np.arange(0,100,0.5)
    #yMC = np.array([f(x).n for x in xMC])
    #yErrMC = np.array([f(x).s for x in xMC])

    print(">>> make plot and combine measurements into bins")
    xx = np.arange(data[xAxis].min()*(1./x_step)//1/(1./x_step),data[xAxis].max()*(1./x_step)//1/(1./x_step), x_step)
    xx = np.append(xx, xx[-1]+x_step)
    xx_centers = []
    yy = []
    yy_err = []
    yy0 = []
    yy0_err = []
    for i in range(0,len(xx)-1):
        dyy = data.loc[(data[xAxis] < xx[i+1]) & (data[xAxis] >= xx[i])][[xAxis, 'y', 'y_Err', 'y0', 'y0_Err']]

        if len(dyy)==0:
            continue

        u_y = np.array([unc.ufloat(y,y_err) for y,y_err in dyy[['y','y_Err']].values])
        u_y0 = np.array([unc.ufloat(y,y_err) for y,y_err in dyy[['y0','y0_Err']].values])

        # 1) simple mean
        yy_avg = u_y.mean()
        yy0_avg = u_y0.mean()

        # 2) or weighted average
        yy_w = np.array([1./(y.s)**2 for y in u_y])
        yy_avg = sum(u_y * yy_w) / sum(yy_w)
        yy0_w = np.array([1./(y.s)**2 for y in u_y0])
        yy0_avg = sum(u_y0 * yy0_w) / sum(yy0_w)

        yy_avg_err = yy_avg.s
        yy_avg = yy_avg.n
        yy0_avg_err = yy0_avg.s
        yy0_avg = yy0_avg.n
        
        if yy_avg == 0 or yy_avg_err/yy_avg > 0.05:
            continue

        xx_centers.append((xx[i] + (xx[i+1] - xx[i]) / 2.))
        yy.append(yy_avg)
        yy_err.append(yy_avg_err)
        yy0.append(yy0_avg)
        yy0_err.append(yy0_avg_err)

    xx = np.array(xx_centers)
    yy = np.array(yy)
    yy_err = np.array(yy_err)
    yy0 = np.array(yy0)
    yy0_err = np.array(yy0_err)

    print(">>> make fit")
    func = linear #pol2 #quad
    popt, pcov = curve_fit(func, xx, yy, sigma=yy_err, absolute_sigma=True)

    perr = np.sqrt(np.diag(pcov))
    params = unc.correlated_values(popt, pcov)

    f = lambda x: func(x, *params)
    xMC = np.arange(0,100,0.5)
    yMC = np.array([f(x).n for x in xMC])
    yErrMC = np.array([f(x).s for x in xMC])

    print(">>> make plots")

    xx_err = np.ones(len(xx_centers))*x_step/2
    
    plt.clf()
    fig = plt.figure(figsize=(10.0,4.0))
    ax1 = fig.add_subplot(111)
    fig.subplots_adjust(left=0.1, right=0.99, top=0.99, bottom=0.16)

    ax1.set_xlabel(xTitle)
    ax1.set_ylabel(yLabel)
    ax1.text(0.74, 0.97, "\\bf{CMS}", verticalalignment='top', transform=ax1.transAxes, weight="bold")
    ax1.text(0.81, 0.97, "\\emph{Work in progress}", verticalalignment='top', transform=ax1.transAxes,style='italic')    
    ax1.text(0.74, 0.89, year, verticalalignment='top', transform=ax1.transAxes,style='italic')    
    
    xMin = min(xx_centers) - x_step/2
    xMax = max(xx_centers) + x_step/2
    xRange = abs(xMin - xMax)

    yMin = 0.95#min(yy - yy_err)
    yMax = 1.05#max(yy + yy_err)
    yRange = abs(yMax - yMin)    

    p4 = ax1.errorbar(xx, yy0, xerr=xx_err, yerr=yy0_err, label="Measurements unc.",
        marker="o", linewidth=0, color="grey", ecolor="grey", elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
        zorder=1)
    
    p3 = ax1.errorbar(xx, yy, xerr=xx_err, yerr=yy_err, label="Measurements",
        marker="o", linewidth=0, color="black", ecolor="black", elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
        zorder=1)
        
    ax1.plot(np.array([xMin-xRange*0.02,xMax+xRange*0.02]),np.array([1.,1.]), "k--", linewidth=1)
    
    p2 = ax1.plot(xMC, yMC, color=colors[0],linestyle="solid",  label="Linear fit", linewidth=1)
    
    p1 = ax1.fill_between(xMC, yMC - yErrMC, yMC + yErrMC,
                     color='grey', alpha=0.2, zorder=1) 
    p1 = ax1.fill(np.NaN, np.NaN, color='grey', alpha=0.2, linewidth=0.)    

    leg_styles = [p4, p3, (p2[0], p1[0])]
    leg_labels = ['Measurements unc.','Measurements', 'Linear fit']
    
    leg = ax1.legend(leg_styles, leg_labels, loc="lower left", ncol=3,
        frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
    leg.get_frame().set_linewidth(0.8)

    ax1.set_ylim([yMin-yRange*0.05, yMax+yRange*0.05])
    ax1.set_xlim([xMin-xRange*0.02, xMax+xRange*0.02])

    outstring = "{0}/{1}_vs_{2}".format(outDir,yAxis,xAxis)
    outstring += year.replace(" ","_")

    plt.savefig(outstring+".png")
    plt.savefig(outstring+".pdf")
    plt.close()

# ------------------------------------------------------------------------------

plot_lowPU = False
if args.xsec:
    print(">>> load low pileup data")
    plot_lowPU = True
    # --- get Z xsec
    data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
    # xsecBB = sum(data_xsec['zDelBB'])/sum(data_xsec['lumiRec'])
    # xsecBE = sum(data_xsec['zDelBE'])/sum(data_xsec['lumiRec'])
    # xsecEE = sum(data_xsec['zDelEE'])/sum(data_xsec['lumiRec'])
    # xsecI = sum(data_xsec['zDelI'])/sum(data_xsec['lumiRec'])

    data_xsec['zDelI'] = data_xsec['zDel'].apply(lambda x: unc.ufloat_fromstr(x).n)

    apply_muon_prefire(data_xsec)
    apply_ECAL_prefire(data_xsec)
    
    xsec = sum(data_xsec['zDelI'])/sum(data_xsec['lumiRec'])

    

# ------------------------------------------------------------------------------
print(">>> load csv file in dataframe")

data_rates = pd.read_csv(args.rates, sep=',')

invalid_runs = {
    275657, 275658, 275659, # Outliers in all those runs of 2016. HFOC was used -> problem there?
    278017, 278018          # More outliers, not clear from where
}
for run in invalid_runs:
    data_rates = data_rates.loc[data_rates['run'] != run]

data_rates['lumiRec'] = data_rates['lumiRec'] * 1000    # convert into /nb

# convert to uncertainties
# data_rates['HLTeffBB'] = data_rates['HLTeffBB'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['HLTeffBE'] = data_rates['HLTeffBE'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['HLTeffEE'] = data_rates['HLTeffEE'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['HLTeffI'] = data_rates['HLTeff'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['SeleffB'] = data_rates['SeleffB'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['SeleffE'] = data_rates['SeleffE'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['SeleffI'] = data_rates['Seleff'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['TrkeffB'] = data_rates['TrkeffB'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['TrkeffE'] = data_rates['TrkeffE'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['TrkeffI'] = data_rates['Trkeff'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['StaeffB'] = data_rates['StaeffB'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['StaeffE'] = data_rates['StaeffE'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['StaeffI'] = data_rates['Staeff'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['GloeffB'] = data_rates['GloeffB'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['GloeffE'] = data_rates['GloeffE'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['GloeffI'] = data_rates['GloeffI'].apply(lambda x: unc.ufloat_fromstr(x))

# data_rates['ZBBeff'] = data_rates['ZBBeff'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['ZBEeff'] = data_rates['ZBEeff'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['ZEEeff'] = data_rates['ZEEeff'].apply(lambda x: unc.ufloat_fromstr(x))
# data_rates['ZIeff'] = data_rates['ZIeff'].apply(lambda x: unc.ufloat_fromstr(x))
# 
# # take uncertainties from uncorrected zDel
# data_rates['zDelBB_mc'] = data_rates['zDelBB_mc'].apply(lambda x: unc.ufloat(x, np.sqrt(x)))
# data_rates['zDelBE_mc'] = data_rates['zDelBE_mc'].apply(lambda x: unc.ufloat(x, np.sqrt(x)))
# data_rates['zDelEE_mc'] = data_rates['zDelEE_mc'].apply(lambda x: unc.ufloat(x, np.sqrt(x)))
# data_rates['zDelI_mc'] = data_rates['zDelI_mc'].apply(lambda x: unc.ufloat(x, np.sqrt(x)))
# 
# data_rates['zDelBB'] = data_rates['zDelBB'].apply(lambda x: unc.ufloat(x, np.sqrt(x)))
# data_rates['zDelBE'] = data_rates['zDelBE'].apply(lambda x: unc.ufloat(x, np.sqrt(x)))
# data_rates['zDelEE'] = data_rates['zDelEE'].apply(lambda x: unc.ufloat(x, np.sqrt(x)))
# data_rates['zDelI'] = data_rates['zDel'].apply(lambda x: unc.ufloat(x, np.sqrt(x)))

data_rates['zDelI'] = data_rates['zDel'].apply(lambda x: unc.ufloat_fromstr(x))

apply_muon_prefire(data_rates)
apply_ECAL_prefire(data_rates)
# apply_pileup_correction(data_rates)
# data_rates['xsecBB_mc'] = data_rates['zDelBB_mc'] / data_rates['lumiRec']
# data_rates['xsecBE_mc'] = data_rates['zDelBE_mc'] / data_rates['lumiRec']
# data_rates['xsecEE_mc'] = data_rates['zDelEE_mc'] / data_rates['lumiRec']
# data_rates['xsecI_mc'] = data_rates['zDelI_mc'] / data_rates['lumiRec']
# data_rates['xsec_mc'] = data_rates['xsecBB_mc'] + data_rates['xsecBE_mc'] + data_rates['xsecEE_mc']

# data_rates['xsecBB'] = data_rates['zDelBB'] / data_rates['lumiRec']
# data_rates['xsecBE'] = data_rates['zDelBE'] / data_rates['lumiRec']
# data_rates['xsecEE'] = data_rates['zDelEE'] / data_rates['lumiRec']
data_rates['xsecI_mc'] = data_rates['zDelI'] / data_rates['lumiRec']
# data_rates['xsec'] = data_rates['xsecBB'] + data_rates['xsecBE'] + data_rates['xsecEE']

data_rates['xsecI'] = data_rates['zDelI'] / data_rates['cIO']**2 / data_rates['lumiRec']

for yy, ylabel, region, mcRes in (
    # ("xsec", "sigma", "", ""),
    # ("xsecBB", "sigma", "BB", ""),
    # ("xsecBE", "sigma", "BE", ""),
    # ("xsecEE", "sigma", "EE", ""),
    ("xsecI_mc", "sigma", "I", ""),
    # # ("xsec_mc", "sigma", "", ""),
    # ("xsecBB_mc", "sigma", "BB", ""),
    # ("xsecBE_mc", "sigma", "BE", ""),
    # ("xsecEE_mc", "sigma", "EE", ""),
    # ("xsecI_mc", "sigma", "I", ""),
    # ("ZBBeff", "Z efficiency","BB", "effBB"),
    # ("ZBEeff", "Z efficiency","BE", "effBE"),
    # ("ZEEeff", "Z efficiency","EE", "effEE"),
    # ("ZIeff", "Z efficiency", "I", "eff"),
    # ('HLTeffBB' ,'Muon HLT efficiency', "BB", "HLTB"),
    # ('HLTeffBE' ,'Muon HLT efficiency', "BE", "HLTE"),
    # ('HLTeffEE' ,'Muon HLT efficiency', "EE", "HLTE"),
    # ('HLTeffI' ,'Muon HLT efficiency', "I", "HLT"),
    # ('SeleffB' ,'Muon selelction efficiency', "B", "SelB"),
    # ('SeleffE' ,'Muon selelction efficiency', "E", "SelE"),
    # ('SeleffI' ,'Muon selelction efficiency', "I", "Sel"),
    # ('GloeffB' ,'Global muon efficiency', "B", "GloB"),
    # ('GloeffE' ,'Global muon efficiency', "E", "GloE"),
    # ('TrkeffB' ,'Muon inner track efficiency', "B", "TrkB"),
    # ('TrkeffE' ,'Muon inner track efficiency', "E", "TrkE"),
    # ('TrkeffI' ,'Muon inner track efficiency', "I", "Trk"),
    # ('StaeffB' ,'Muon standalone efficiency', "B", "StaB"),
    # ('StaeffE' ,'Muon standalone efficiency', "E", "StaE"),
    # ('StaeffI' ,'Muon standalone efficiency', "I", "Sta"),
):
    # # # single eras
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 B", run_range=(272007,275376), normalized=False)
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 C", run_range=(275657,276283), normalized=False)
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 D", run_range=(276315,276811), normalized=False)
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 E", run_range=(276831,277420), normalized=False)
    # # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 F", run_range=(277772,278808), normalized=False)
    # # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 G", run_range=(278820,280385), normalized=False)
    # # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 H", run_range=(280919,284044), normalized=False)
    # # 
    # # # 2017
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 B", title="corrected", run_range=(297046,299329), normalized=False)
    # # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 C", title="corrected", run_range=(299368,302029), normalized=False)
    # # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 D", title="corrected", run_range=(302030,303434), normalized=False)
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 E", title="corrected", run_range=(303434,304797), normalized=False)
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 F", title="corrected", run_range=(305040,306462), normalized=False)
    # # 
    # # # ## 2018
    # # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 A", run_range=(315252,316995), normalized=False)
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 B", run_range=(317080,319310), normalized=False)
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 C", run_range=(319337,320065), normalized=False)
    # # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 D", run_range=(320673,325175), normalized=False)

    # # # total 2016
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016", run_range=(272007,294645), normalized=True)
    # # total 2016 pre VFP
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 pre VFP", run_range=(272007,278769), normalized=True)
    # total 2016 post VFP
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 post VFP", run_range=(278769,294645), normalized=True)
    # # total 2017
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2017", run_range=(297046,306462), normalized=True)
    # # ## total 2018
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018", run_range=(315252,325175), normalized=True)
