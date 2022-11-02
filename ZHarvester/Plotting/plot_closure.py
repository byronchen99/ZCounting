import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

import json
import numpy as np
import uncertainties as unc
import pdb
from scipy.optimize import curve_fit
import pandas as pd
import os, sys
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))

from ZUtils.python.utils import linear

textsize = 16
markersize = 4
ylabel = "Correlation factor"

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

colors = ["#e74c3c","#2980b9","#27ae60","#f1c40f","#8e44ad"]

xLabels = {
    "nPV": "Number of primary vertices",
    "nPU": "Number of pileup events",
    "pT": "Muon pT",
    "absEta1": "Muon $|\\eta|$",
    "absEta2": "Muon $|\\eta|$",
    "eta": "Muon eta",
    "dz": "$\\frac{z(\\mu^+)+z(\\mu^-)}{2}$ [cm]"
}
triggers = {
    "2016preVFP": "HLT\_IsoMu24 or HLT\_IsoTkMu24",
    "2016postVFP": "HLT\_IsoMu24 or HLT\_IsoTkMu24",
    "2017": "HLT\_IsoMu24",
    "2018": "HLT\_IsoMu24",
    "Summer16preVFP": "HLT\_IsoMu24 or HLT\_IsoTkMu24",
    "Summer16postVFP": "HLT\_IsoMu24 or HLT\_IsoTkMu24",
    "Fall17": "HLT\_IsoMu24",
    "Autumn18": "HLT\_IsoMu24"
}
yearnames = {
    "2016preVFP": "2016 preVFP",
    "2016postVFP": "2016 postVFP",
    "2017": "2017",
    "2018": "2018",
    "Summer16preVFP": "2016 preVFP",
    "Summer16postVFP": "2016 postVFP",
    "Fall17": "2017",
    "Autumn18": "2018"
}
xVar = "nPU"

massMin_ = 60
massMax_ = 120

ptCut_ = 25.
delR_ = 0.0

massMinSta_ = 50
massMaxSta_ = 130

ptCutSta_= 25

nBinsX = 18
minX = -0.5
maxX = 71.5
widthX = (maxX - minX )/nBinsX
useCorrIO=True

def preparePlot(yLabel, observable=""):

    xLabel = xLabels.get(observable,"")

    plt.clf()
    fig = plt.figure()
    fig.subplots_adjust(left=0.15, right=0.98, top=0.97, bottom=0.125)
    ax = fig.add_subplot(111)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    return ax

def finalizePlot(ax, legend="best", prefix="", cmsTag="", year="", nameTag="", region="", yRange=None):
    # ax.set_xlim([1,74])

    # xT = 0.02
    # yT = 0.22

    xT = 0.02
    yT = 0.98

    reg = ""
    if region =="I":
        reg = "Inclusive"
    if region =="BB":
        reg = "Barrel-Barrel"
    if region =="BE":
        reg = "Barrel-Endcap"
    if region =="EE":
        reg = "Endcap-Endcap"

    #if "HLT" in prefix:
    #    xT = 0.42
    #    # ax.text(0.02, 0.02, triggers[year[1]], verticalalignment='bottom', transform=ax.transAxes)
    #else:
    #    xT = 0.30

    ax.text(xT, yT, "\\bf{CMS}", verticalalignment='top', transform=ax.transAxes, weight="bold")
    ax.text(xT+0.11, yT, "\\emph{Simulation}", verticalalignment='top', transform=ax.transAxes,style='italic')
    ax.text(xT, yT-0.07, "\\emph{Work in progress}", verticalalignment='top', transform=ax.transAxes,style='italic')    
    ax.text(0.98, yT, yearnames[year[1]], verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)


    leg = ax.legend(loc=legend, frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
    leg.get_frame().set_linewidth(0.8)

    ax.set_xlim(rangeX)
    # if year == "2016preVFP":
    #     ax.set_ylim(0.3,1.1)
    # elif year == "2016postVFP":
    #     ax.set_ylim(0.84,1.01)
    # elif year == "2017":
    #     ax.set_ylim(0.96,1.01)
    # elif year == "2018":
    if yRange:
        ax.set_ylim(yRange)

    outputname = "MCClosure_{0}_{1}_{2}_{3}".format(prefix, year[1], nameTag, region)
    print(outputname)
    plt.savefig(outputname+".png")
    plt.savefig(outputname+".pdf")
    plt.close()


for year in [
   ["2016preVFP","Summer16preVFP"],
   ["2016postVFP","Summer16postVFP"],
   ["2017","Fall17"],
   ["2018","Autumn18"]
    ]:
        
    if year[0] == "2016preVFP":
        rangeX = [0,54]
        rangeY = [0.9971,1.0071]

    elif year[0] == "2016postVFP":
        rangeX = [0,60]
        rangeY = [0.9971,1.0071]
    else:
        rangeX = [0,60]
        rangeY = [0.9971,1.0071]

    rangeX = [0,72]

    yearname = yearnames[year[0]]

    trigger = triggers[year[0]]

    xlabel = xLabels[xVar]

    # ----------------------------------------------------------------------
    # Get correlation factor on outer-inner track
    if useCorrIO:
        # tfileIO = ROOT.TFile("/nfs/dust/cms/user/dwalter/CMSSW_10_2_20_UL/src/potato-zcount/zcount/output/ZMCClosureGen-V17_01-DYFlat-d20220331-t174316/HistogramerMCClosureGen-V17_01-{0}-DYJetsToLL_M_50_LO.root".format(year[1]),"READ")

        tfileIO = ROOT.TFile("/nfs/dust/cms/user/dwalter/ZCounting/CMSSW_12_4_9/src/ZCounting/ZHarvester/res/mcCorrections/correlations_V17_30/c_nPU_{0}.root".format(year[0]),"READ")
                        
        def hist2array(name):
            _h = tfileIO.Get(name)

            _a = np.array([_h.GetBinContent(x) for x in range(1, _h.GetNbinsX()+1)])
            _b = np.array([sum(_a[i*int(widthX):(i+1)*int(widthX)])/int(widthX) for i in range(int(len(_a)/widthX))])
            return _b

        arr_cIO = hist2array("cIO_I")
        arr_cID = hist2array("cID_I")

    # tfile = ROOT.TFile("/afs/desy.de/user/d/dwalter/nfsHome/data/Lumi/V17/TTrees/ZCountingAll-V17_01-{0}-DYJetsToLL_M_50_LO.root".format(year[1]),"READ")
    tfile = ROOT.TFile("/afs/desy.de/user/d/dwalter/nfsHome/data/Lumi/V17_30/TTrees/ZCountingAll-V17_30-{0}-DYJetsToLL_M_50_LO.root".format(year[1]),"READ")
    # tfile = ROOT.TFile("/nfs/dust/cms/user/dwalter/CMSSW_10_2_20_UL/src/potato-zcount/zcount/output/ZCountingAll-V17_28-DYFlatV17-d20221031-t141632/ZCountingAll-V17_28-{0}-DYJetsToLL_M_50_LO.root".format(year[1]),"READ")

    print("load events in arrays")
   
    tTruth = tfile.Get("Truth")        
    tHLT = tfile.Get("HLT")
    tSel = tfile.Get("Sel")
    tTrk = tfile.Get("Trk")
    
    hTruth = ROOT.TH1D("htruth","Truth ", nBinsX,minX,maxX)
    hlt0 = ROOT.TH1D("hlt0","HLT 0", nBinsX,minX,maxX)
    hlt1 = ROOT.TH1D("hlt1","HLT 1", nBinsX,minX,maxX)
    hlt2 = ROOT.TH1D("hlt2","HLT 2", nBinsX,minX,maxX)
    selFail = ROOT.TH1D("selFail","SelFail", nBinsX,minX,maxX)
    
    trkPass = ROOT.TH1D("trkPass","TrkPass", nBinsX,minX,maxX)
    trkFail = ROOT.TH1D("trkFail","TrkFail", nBinsX,minX,maxX)    

    # # for higher order corrections on the efficiency
    # trkPass2 = ROOT.TH1D("trkPass2","TrkPass2", nBinsX,minX,maxX)
    # trkPass1 = ROOT.TH1D("trkPass1","TrkPass1", nBinsX,minX,maxX)

    tTruth.Draw("{0}>>htruth".format(xVar))

    tHLT.Draw("{0}>>hlt0".format(xVar),"match1 && match2 && pass==0    && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMin_,massMax_, ptCut_, delR_))
    tHLT.Draw("{0}>>hlt1".format(xVar),"match1 && match2 && pass==1    && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMin_,massMax_, ptCut_, delR_))
    tHLT.Draw("{0}>>hlt2".format(xVar),"match1 && match2 && pass==2    && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMin_,massMax_, ptCut_, delR_))
    tSel.Draw("{0}>>selFail".format(xVar),"match1 && match2 && pass==0 && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMin_,massMax_, ptCut_, delR_))

    tTrk.Draw("{0}>>trkFail".format(xVar),"match1 && match2 && pass==0 && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMinSta_,massMaxSta_, ptCutSta_, delR_))
    # tTrk.Draw("{0}>>trkPass".format(xVar),"match1 && match2 && pass>=1 && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMinSta_,massMaxSta_, ptCutSta_, delR_))

    tTrk.Draw("{0}>>trkPass".format(xVar),"match1 && match2 && pass>=1 && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMinSta_,massMaxSta_, ptCutSta_, delR_))
    #  \
    #     && massTrk>{4} && massTrk<{5} && ptTrk>{6}".format(massMinSta_,massMaxSta_, ptCutSta_, delR_, massMin_,massMax_, ptCutTrk_))

    # tTrk.Draw("{0}>>trkPass2".format(xVar),"match1 && match2 && pass==2 && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMinSta_,massMaxSta_, ptCut_, delR_))
    # tTrk.Draw("{0}>>trkPass1".format(xVar),"match1 && match2 && pass==1 && mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2} && delR > {3}".format(massMinSta_,massMaxSta_, ptCut_, delR_))

    def hist2array(_h):
        return np.array([ unc.ufloat(_h.GetBinContent(x), np.sqrt(_h.GetBinContent(x))) for x in range(1, _h.GetNbinsX()+1)])

    hists = [
        hTruth,
        hlt0, hlt1, hlt2,
        selFail,
        trkFail, 
        trkPass,
        # trkPass2,
        # trkPass1
    ]
    keys = [
        "Truth",
        "HLT0", "HLT1", "HLT2",
        "SelFail",
        "TrkFail", 
        "TrkPass",
        # "TrkPass2",
        # "TrkPass1"         
    ]
    arrs = [hist2array(key) for key in hists]    

    keys.append("x")
    # arrs.append( (np.linspace(minX,maxX,nBinsX+widthX)+widthX/2.)[:-1] )
    arrs.append(np.arange(minX+0.5,maxX+0.5,widthX))#+widthX/2.)

 
    print("make dataframe")

    df = pd.DataFrame(np.array(arrs).T, columns=keys)

    # Di-muon events

    # calculate efficiencies
    # tnp efficiencies
    df["effHLT"] = (2*df["HLT2"]) / (2*df["HLT2"] + df["HLT1"])
    df["effSel"] = (2*df["HLT2"] + df["HLT1"]) / (2*df["HLT2"] + df["HLT1"] + df["SelFail"])
    
    df["effTrk"] = df["TrkPass"] / (df["TrkPass"] + df["TrkFail"])
    # df["effTrk"] = (df["TrkPass1"] + df["TrkPass2"] - (df["TrkPass2"]*df["TrkFail"])/df["TrkPass1"] ) / (df["TrkPass1"] + df["TrkPass2"] + df["TrkFail"])
    
    df["effID"] = df["effSel"] * df["effTrk"]

    df["cHLT"] = 4 * (df["HLT2"]+df["HLT1"]+df["HLT0"]) * df["HLT2"] / (df["HLT1"] + 2*df["HLT2"])**2
            
    # ----- closure test for Z -> mumu
    # calculate true numbers at each step

    # from HLTTruth to RecoTruth
    df["NZhlt_reco"] = (df["HLT2"]+df["HLT1"]/2.)**2 / df["HLT2"]
    
    # from HLTTruth to GloTruth
    df["NZhlt_glo"] = df["NZhlt_reco"] / df["effSel"]**2 
    
    # from HLT to Truth
    df["NZhlt"] = df["NZhlt_reco"] / (df["effID"])**2
    
    df["a"] = df["Truth"]        
    # apply acceptance corrections
            
    # di muon acceptance
    acceptanceFactor = (df["NZhlt"]).sum() / df["Truth"].sum()
    print("Acceptance factor (w/o correlation factors) {0}".format(acceptanceFactor))

    acceptanceFactor = (df["NZhlt"] * df["cHLT"]).sum() / df["Truth"].sum()
    print("Acceptance factor for (w/ c_HLT ) {0}".format(acceptanceFactor))

    if useCorrIO:
        df["cID"] = arr_cIO[:nBinsX]
        acceptanceFactor = (df["NZhlt"] * df["cHLT"] * df["cID"]).sum() / df["Truth"].sum()
        print("Acceptance factor (w/ c_HLT c_ID) {0}".format(acceptanceFactor))

        df["cIO"] = arr_cIO[:nBinsX]
        acceptanceFactor = (df["NZhlt"] * df["cHLT"] * df["cID"] * df["cIO"]**2).sum() / df["Truth"].sum()
        print("Acceptance factor (w/ c_io^2 c_ID c_HLT) {0}".format(acceptanceFactor))



    
    # acceptanceFactor = 0.306180735128

    df["RZhlt"] = df["NZhlt"] / (df["Truth"]*acceptanceFactor)
    
    xx = df["x"].values
    
    # if year[1] == "Summer16preVFP":
    #     rangeY=[0.82,1.025]
    # elif year[1] == "Summer16postVFP":
    #     rangeY=[0.95,1.025]
    # elif year[1] == "Fall17":
    #     rangeY = [0.98, 1.01]
    # else:
    rangeY=[0.97,1.0175]
    
    ax = preparePlot("Measurement / Truth", xVar)
    ax.plot([min(xx), max(xx+widthX)], [1,1], "k--")

    xx0 = np.append(xx, xx[-1]+widthX)
    yy0 = df["RZhlt"].apply(lambda x: x.n).values
    yy0 = np.append(yy0, yy0[-1])

    ax.plot(xx0, yy0, 
        label="$N_\\mathrm{Z}$", drawstyle='steps-post', color=colors[1], zorder=0)                  
    # ax.plot(xx, df["RZhlt_corr1"].values, label="$N_\\mathrm{Z}\ C_\\mathrm{HLT}$", drawstyle='steps', color=colors[2])    

    # do_fits = False
    # if do_fits:
    #     # perform fit to quantify linearity
    #     func = linear
    #     popt, pcov = curve_fit(func, xx[1:], df["RZhlt_corr1"].apply(lambda x: x.n).values[1:])

    #     perr = np.sqrt(np.diag(pcov))
    #     # correlated values
    #     params = unc.correlated_values(popt, pcov)

    #     yy_corr1 = np.array([func(x, *params).n for x in xx])

    #     func = linear
    #     popt, pcov = curve_fit(func, xx[1:], df["RZhlt"].apply(lambda x: x.n).values[1:])

    #     perr = np.sqrt(np.diag(pcov))
    #     # correlated values
    #     params = unc.correlated_values(popt, pcov)

    #     yy = np.array([func(x, *params).n for x in xx])

    #     ax.plot(xx, yy, color=colors[1])    
    #     ax.plot(xx, yy_corr1, color=colors[2])      

    if useCorrIO:
        
        df["RZhlt_corrIO"] = df["RZhlt"] * df["cIO"]**2        
        df["RZhlt_corrIOandID"] = df["RZhlt_corrIO"] * df["cID"]
        df["RZhlt_corrIOandIDandHLT"] = df["RZhlt_corrIOandID"] * df["cHLT"]

        xx0 = np.append(xx, xx[-1]+widthX)
        yy0 = df["RZhlt_corrIO"].apply(lambda x: x.n).values
        yy0 = np.append(yy0, yy0[-1])

        ax.plot(xx0, yy0, 
            label="$N_\\mathrm{Z}\ C_\\mathrm{IO}^2$", drawstyle='steps-post', 
            color=colors[2], zorder=1)    
        
 
        yy0 = df["RZhlt_corrIOandID"].apply(lambda x: x.n).values
        yy0 = np.append(yy0, yy0[-1])

        ax.plot(xx0, yy0, 
            label="$N_\\mathrm{Z}\ C_\\mathrm{IO}^2 C_\\mathrm{ID}$", 
            drawstyle='steps-post', 
            color=colors[3], zorder=1)            

        xx0 = xx+widthX/2.
        yy0 = df["RZhlt_corrIOandIDandHLT"].apply(lambda x: x.n).values
        yyerr0 = df["RZhlt_corrIOandIDandHLT"].apply(lambda x: x.s).values
        xxerr0 = np.ones(len(xx))*widthX/2.

        ax.errorbar(xx0, yy0, yerr=yyerr0, xerr=xxerr0,
            label="$N_\\mathrm{Z}\ C_\\mathrm{IO}^2 C_\\mathrm{ID} C_\\mathrm{HLT}$", 
            #drawstyle='steps-post', 
            color=colors[0], linestyle="",marker=".", zorder=2)            

        # --> Do fit
        print(">>> make fit")
        func = linear #pol2 #quad
        popt, pcov = curve_fit(func, xx0, yy0, sigma=yyerr0, absolute_sigma=True)

        perr = np.sqrt(np.diag(pcov))
        params = unc.correlated_values(popt, pcov)

        f = lambda x: func(x, *params)
        xFit = np.arange(0,100,0.5)
        yFit = np.array([f(x).n for x in xFit])
        yErrFit = np.array([f(x).s for x in xFit])

        ax.plot(xFit, yFit, color=colors[0])   
        # <--

    else:
        df["RZhlt_corrHLT"] = df["RZhlt"] * df["cHLT"]

        ax.errorbar(xx+widthX/2., df["RZhlt_corrHLT"].apply(lambda x: x.n).values, 
            yerr=df["RZhlt_corrHLT"].apply(lambda x: x.s).values, xerr=np.ones(len(xx))*widthX/2.,
            label="$N_\\mathrm{Z}\ C_\\mathrm{HLT}$", 
            #drawstyle='steps-post', 
            color=colors[0], linestyle="",marker=".", zorder=2)              

    finalizePlot(ax, "lower left", "HLT", region="I", cmsTag="center left", year=year, nameTag="", yRange=rangeY)                
        
    # for c, yRangeEff in (
    #     ("HLT", [0.70,1.03]),
    #     ("ID",  [0.90,1.02]),
    #     ("Sel", [0.90,1.01]),
    #     ("Glo", [0.95,1.005]),
    #     ("Sta", [0.95,1.005]),
    #     ("Trk", [0.7,1.005]),
    # ):
    #     if "eff"+c not in df.keys():
    #         continue;
    #     ax = preparePlot("Efficiency", xVar)
    #     ax.plot([min(xx), max(xx)], [1,1], "k--")
    # 
    #     if dfInfo is not None and "eff{0}".format(c) in dfInfo.keys():
    #         ax.errorbar(dfInfo["x"].values, dfInfo["eff{0}".format(c)].values, yerr=dfInfo["eff{0}_err".format(c)].values, xerr=dfInfo["xErr"].values, 
    #             color="black", linestyle="",marker=".")       
    # 
    #     ax.plot(xx, df["eff{}".format(c)].values, "k-", label="$\epsilon_\mathrm{"+c+"}$ (I)", drawstyle='steps')
    #

    # for c, yRangeC, doFit in (
    #     ("HLT", [0.995,1.01], True),
    # ):
    #     ax = preparePlot("Correlation coefficient", xVar)
    #     ax.plot([min(xx), max(xx)], [1,1], "k--")
    #     ax.step(xx, df["c{}".format(c)].values, "k-", label="Simulation", where="mid")
    # 
    #     if dfInfoCorrHLT is not None and "c" in dfInfoCorrHLT.keys():
    #         ax.errorbar(dfInfoCorrHLT["x"].values, dfInfoCorrHLT["c".format(c)].values, 
    #             yerr=dfInfoCorrHLT["c_err".format(c)].values, 
    #             xerr=(dfInfoCorrHLT["xLo".format(c)].values,dfInfoCorrHLT["xHi".format(c)].values),
    #             color="black", linestyle="",marker=".", label="Data")
    # 
    #     if doFit:
    #         func = linear
    #         popt, pcov = curve_fit(func, xx, df["c{}".format(c)].values)
    # 
    #         perr = np.sqrt(np.diag(pcov))
    #         # correlated values
    #         params = unc.correlated_values(popt, pcov)
    # 
    #         yy = np.array([func(x, *params).n for x in xx])
    #         yyErr = np.array([func(x, *params).s for x in xx])
    # 
    #         ax.plot(xx, yy, color="black")    
    #         ax.fill_between(xx, yy - yyErr, yy + yyErr, 
    #             # step = 'pre',
    #             color='grey', alpha=0.5, zorder=0)   
    # 
    #     finalizePlot(ax, "upper left", "c{}".format(c),region="All", year=year, nameTag="", yRange=yRangeC, simulation=False)    
