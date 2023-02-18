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
import matplotlib.pyplot as plt
import argparse

sys.path.append(os.getcwd())
os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))

from ZUtils.python.utils import linear

from common import parsing, plotting
from common.logging import child_logger
log = child_logger(__name__)

parser = parsing.parser_plot()

args = parser.parse_args()

output = args.output
if not os.path.isdir(output):
    os.mkdir(output)

colors, textsize, labelsize, markersize = plotting.set_matplotlib_style()

ylabel = "Correlation factor"

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
    "2022": "HLT\_IsoMu24",
    "Summer16preVFP": "HLT\_IsoMu24 or HLT\_IsoTkMu24",
    "Summer16postVFP": "HLT\_IsoMu24 or HLT\_IsoTkMu24",
    "Fall17": "HLT\_IsoMu24",
    "Autumn18": "HLT\_IsoMu24",
    "Winter22": "HLT\_IsoMu24"
}
yearnames = {
    "2016preVFP": "2016 preVFP",
    "2016postVFP": "2016 postVFP",
    "2017": "2017",
    "2018": "2018",
    "2022": "2022",
    "Summer16preVFP": "2016 preVFP",
    "Summer16postVFP": "2016 postVFP",
    "Fall17": "2017",
    "Autumn18": "2018",
    "Winter22": "2022"
}

ptCut = 27.
etaCut = 2.4

delR_ = 0.0

massMin_ = 66
massMax_ = 116

massMinSta_ = 50
massMaxSta_ = 130

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

def finalizePlot(ax, legend="best", prefix="", cmsTag="", year="", nameTag="", region="", yRange=None, ncol=2):
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
    # ax.text(xT, yT-0.07, "\\emph{Preliminary}", verticalalignment='top', transform=ax.transAxes,style='italic')    
    ax.text(0.98, yT, yearnames[year[1]], verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)


    leg = ax.legend(loc=legend, frameon=True, ncol=ncol,framealpha=1.0, fancybox=False, edgecolor="black")
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
    print(output+"/"+outputname+".png")
    plt.savefig(output+"/"+outputname+".png")
    plt.savefig(output+"/"+outputname+".pdf")
    plt.close()

for xVar in [
    "nPU", 
    "nPV",
]:
    for year in [
        # ["2016preVFP","Summer16preVFP"],
        # ["2016postVFP","Summer16postVFP"],
        # ["2017","Fall17"],
        # ["2018","Autumn18"],
        ["2022","Winter22"],
    ]:

        if year[0] == "2022":
            nBinsX = 70
            minX = -0.5
            maxX = 69.5
            ptCut = 27
        else:
            nBinsX = 75
            minX = -0.5
            maxX = 74.5
            ptCut = 25

        widthX = (maxX - minX )/nBinsX
            
        if year[0] == "2016preVFP":
            rangeX = [0,54]
            rangeY = [0.9971,1.0071]

        elif year[0] == "2016postVFP":
            rangeX = [0,60]
            rangeY = [0.9971,1.0071]
        else:
            rangeX = [0,60]
            
        rangeY = [0.991,1.029]
        rangeX = [0,72]

        yearname = yearnames[year[0]]

        trigger = triggers[year[0]]

        xlabel = xLabels[xVar]

        # ----------------------------------------------------------------------
        # Get correlation factor on outer-inner track
        if useCorrIO:
            # tfileIO = ROOT.TFile("/nfs/dust/cms/user/dwalter/ZCounting/CMSSW_12_4_9/src/ZCounting/ZHarvester/res/mcCorrections/correlations_V17_38/c_{0}_{1}.root".format(xVar, year[0]),"READ")

            tfileIO = ROOT.TFile("/eos/home-d/dwalter/Lumi/ZCounting_potato/CMSSW_10_2_20_UL/src/potato-zcount/plots/MCClosureGen-V19_07-d20230218-t160908/c_{0}_{1}.root".format(xVar, year[0]),"READ")

            def hist2array(name):
                _h = tfileIO.Get(name)

                _a = np.array([_h.GetBinContent(x) for x in range(1, _h.GetNbinsX()+1)])
                _b = np.array([sum(_a[i*int(widthX):(i+1)*int(widthX)])/int(widthX) for i in range(int(len(_a)/widthX))])
                return _b

            arr_cIO = hist2array("cIO_I")
            arr_CID = hist2array("CID_I")

        # tfile = ROOT.TFile("/eos/home-d/dwalter/Lumi/ZCounting_potato/CMSSW_10_2_20_UL/src/potato-zcount/None/ZCountingInOut-V19_01-{0}-DYJetsToLL_M_50_LO-jobI-d20230214-t190442.root".format(year[1]),"READ")
        # tfile = ROOT.TFile("/eos/home-d/dwalter/Lumi/ZCounting_potato/CMSSW_10_2_20_UL/src/potato-zcount/None/ZCountingInOut-V19_02-{0}-DYJetsToLL_M_50_LO-jobI-d20230214-t193941.root".format(year[1]),"READ")
        tfile = ROOT.TFile("/eos/cms/store/group/comm_luminosity/ZCounting/2022/SignalTemplates/ZCountingInOut-V19_07-{0}-DYJetsToLL_M_50_LO.root".format(year[1]),"READ")

        print("load events in arrays")
    
        tTruth = tfile.Get("Truth")        
        tHLT = tfile.Get("HLT")
        tID = tfile.Get("ID")
        tGlo = tfile.Get("Glo")
        tSta = tfile.Get("Sta")
        
        # gen cuts
        hTruth = ROOT.TH1D("htruth","Truth ", nBinsX,minX,maxX)

        hlt0Gen = ROOT.TH1D("hlt0Gen","HLT 0", nBinsX,minX,maxX)
        hlt1Gen = ROOT.TH1D("hlt1Gen","HLT 1", nBinsX,minX,maxX)
        hlt2Gen = ROOT.TH1D("hlt2Gen","HLT 2", nBinsX,minX,maxX)
        
        idFailGen = ROOT.TH1D("idFailGen","idFail", nBinsX,minX,maxX)
        
        gloPassGen = ROOT.TH1D("gloPassGen","glo Pass gen", nBinsX,minX,maxX)
        gloFailGen = ROOT.TH1D("gloFailGen","glo Fail gen", nBinsX,minX,maxX)    

        staPassGen = ROOT.TH1D("staPassGen","sta Pass gen", nBinsX,minX,maxX)
        staFailGen = ROOT.TH1D("staFailGen","sta Fail gen", nBinsX,minX,maxX)

        genCuts = f"ptGen1 > {ptCut} && ptGen2 > {ptCut} && abs(etaGen1) < {etaCut} && abs(etaGen2) < {etaCut} && massGen > {massMin_} && massGen < {massMax_}"

        tTruth.Draw("{0}>>htruth".format(xVar),genCuts)

        tHLT.Draw("{0}>>hlt0Gen".format(xVar),"match1 && match2 && pass==0    && {0}".format(genCuts))
        tHLT.Draw("{0}>>hlt1Gen".format(xVar),"match1 && match2 && pass==1    && {0}".format(genCuts))
        tHLT.Draw("{0}>>hlt2Gen".format(xVar),"match1 && match2 && pass==2    && {0}".format(genCuts))

        tID.Draw("{0}>>idFailGen".format(xVar),"match1 && match2 && pass==0   && {0}".format(genCuts))

        tGlo.Draw("{0}>>gloFailGen".format(xVar),"match1 && match2 && pass==0 && {0}".format(genCuts))
        tGlo.Draw("{0}>>gloPassGen".format(xVar),"match1 && match2 && pass>=1 && {0}".format(genCuts))

        tSta.Draw("{0}>>staFailGen".format(xVar),"match1 && match2 && pass==0 && {0}".format(genCuts))
        tSta.Draw("{0}>>staPassGen".format(xVar),"match1 && match2 && pass>=1 && {0}".format(genCuts))


        # reco cuts
        hlt0 = ROOT.TH1D("hlt0","HLT 0", nBinsX,minX,maxX)
        hlt1 = ROOT.TH1D("hlt1","HLT 1", nBinsX,minX,maxX)
        hlt2 = ROOT.TH1D("hlt2","HLT 2", nBinsX,minX,maxX)

        idFail = ROOT.TH1D("idFail","idFail", nBinsX,minX,maxX)
        
        gloPass = ROOT.TH1D("gloPass","gloPass", nBinsX,minX,maxX)
        gloFail = ROOT.TH1D("gloFail","gloFail", nBinsX,minX,maxX)    

        staPass = ROOT.TH1D("staPass","staPass", nBinsX,minX,maxX)
        staFail = ROOT.TH1D("staFail","staFail", nBinsX,minX,maxX)

        recoCuts    = "mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2}".format(massMin_,massMax_, ptCut)
        recoCutsSta = "mass>{0} && mass<{1} && ptTag > {2} && ptProbe > {2}".format(massMinSta_,massMaxSta_, ptCut)

        tHLT.Draw("{0}>>hlt0".format(xVar),"match1 && match2 && pass==0    && {0}".format(recoCuts))
        tHLT.Draw("{0}>>hlt1".format(xVar),"match1 && match2 && pass==1    && {0}".format(recoCuts))
        tHLT.Draw("{0}>>hlt2".format(xVar),"match1 && match2 && pass==2    && {0}".format(recoCuts))

        tID.Draw("{0}>>idFail".format(xVar),"match1 && match2 && pass==0   && {0}".format(recoCuts))

        tGlo.Draw("{0}>>gloFail".format(xVar),"match1 && match2 && pass==0 && {0}".format(recoCutsSta))
        tGlo.Draw("{0}>>gloPass".format(xVar),"match1 && match2 && pass>=1 && {0}".format(recoCutsSta))

        tSta.Draw("{0}>>staFail".format(xVar),"match1 && match2 && pass==0 && {0}".format(recoCuts))
        tSta.Draw("{0}>>staPass".format(xVar),"match1 && match2 && pass>=1 && {0}".format(recoCuts))


        def hist2array(_h):
            return np.array([ unc.ufloat(_h.GetBinContent(x), np.sqrt(_h.GetBinContent(x))) for x in range(1, _h.GetNbinsX()+1)])

        hists = [
            hTruth,
            hlt0Gen, hlt1Gen, hlt2Gen,
            idFailGen,
            gloFailGen, 
            gloPassGen,
            staFailGen, 
            staPassGen,
            hlt0, hlt1, hlt2,
            idFail,
            gloFail, 
            gloPass,
            staFail, 
            staPass,
        ]
        keys = [
            "Truth",
            "HLT0Gen", "HLT1Gen", "HLT2Gen",
            "idFailGen",
            "gloFailGen", 
            "gloPassGen", 
            "staFailGen", 
            "staPassGen",       
            "HLT0", "HLT1", "HLT2",
            "idFail",
            "gloFail", 
            "gloPass", 
            "staFail", 
            "staPass",       
        ]
        arrs = [hist2array(key) for key in hists]    

        keys.append("x")
        # arrs.append( (np.linspace(minX,maxX,nBinsX+widthX)+widthX/2.)[:-1] )
        arrs.append(np.arange(minX+0.5,maxX+0.5,widthX))#+widthX/2.)

    
        print("make dataframe")

        df = pd.DataFrame(np.array(arrs).T, columns=keys)

        if useCorrIO:
            df["CIO"] = arr_cIO[:nBinsX]
            df["CID"] = arr_CID[:nBinsX]

        # Di-muon events - using gen cuts
        df["effIDGen"] = (2*df["HLT2Gen"] + df["HLT1Gen"]) / (2*df["HLT2Gen"] + df["HLT1Gen"] + df["idFailGen"])
        df["effGloGen"] = df["gloPassGen"] / (df["gloPassGen"] + df["gloFailGen"])
        df["effStaGen"] = df["staPassGen"] / (df["staPassGen"] + df["staFailGen"])
        
        df["effMuGen"] = df["effIDGen"] * df["effGloGen"] * df["effStaGen"]

        df["CHLTGen"] = 4 * (df["HLT2Gen"]+df["HLT1Gen"]+df["HLT0Gen"]) * df["HLT2Gen"] / (df["HLT1Gen"] + 2*df["HLT2Gen"])**2
    
        df["NZ_recoGen"] = (df["HLT2Gen"]+df["HLT1Gen"]/2.)**2 / df["HLT2Gen"]
        df["NZGen"] = df["NZ_recoGen"] / (df["effMuGen"])**2

        df["RZGen"] = df["NZGen"] / (df["Truth"])

        df.at[0,"RZGen"] = unc.ufloat(1.0,0.0)
        df.at[0,"CHLTGen"] = unc.ufloat(1.0,0.0)

        # ----- make plot for gen cuts
        xx = df["x"].values
        xx0 = np.append(xx, xx[-1]+widthX)

        ax = preparePlot("Measurement / Truth", xVar)
        ax.plot([min(xx), max(xx+widthX)], [1,1], "k--")

        yy0 = df["RZGen"].apply(lambda x: x.n).values
        yy0 = np.append(yy0, yy0[-1])
        ax.plot(xx0, yy0, label="$N^\\mathrm{Z}$", drawstyle='steps-mid', color=colors[1], zorder=0)                  

        df["RZGen_corrIO"] = df["RZGen"] * df["CIO"]**2
        df["RZGen_corrIOandID"] = df["RZGen_corrIO"] * df["CID"]
        df["RZGen_corrIOandIDandHLT"] = df["RZGen_corrIOandID"] * df["CHLTGen"]

        yy0 = df["RZGen_corrIO"].apply(lambda x: x.n).values
        yy0 = np.append(yy0, yy0[-1])

        ax.plot(xx0, yy0, label="$N^\\mathrm{Z}\ C_\\mathrm{IO}^2$", drawstyle='steps-mid', color=colors[2], zorder=1)    
            
        yy0 = df["RZGen_corrIOandID"].apply(lambda x: x.n).values
        yy0 = np.append(yy0, yy0[-1])

        ax.plot(xx0, yy0, label="$N^\\mathrm{Z}\ C_\\mathrm{IO}^2 C_\\mathrm{ID}$", drawstyle='steps-mid', color=colors[3], zorder=2)            

        # xx0 = xx+widthX/2.
        yy0 = df["RZGen_corrIOandIDandHLT"].apply(lambda x: x.n).values
        yy0 = np.append(yy0, yy0[-1])

        ax.plot(xx0, yy0, label="$N^\\mathrm{Z}\ C_\\mathrm{IO}^2 C_\\mathrm{ID} C_\\mathrm{HLT}$", drawstyle='steps-mid', color=colors[0], zorder=3)   

        finalizePlot(ax, "lower left", "Gen", region="I", cmsTag="center left", year=year, nameTag=xVar, yRange=[0.971,1.009])                


        # calculate efficiencies with reco cuts
        # tnp efficiencies
        df["effID"] = (2*df["HLT2"] + df["HLT1"]) / (2*df["HLT2"] + df["HLT1"] + df["idFail"])
        df["effGlo"] = df["gloPass"] / (df["gloPass"] + df["gloFail"])
        df["effSta"] = df["staPass"] / (df["staPass"] + df["staFail"])
        
        df["effMu"] = df["effID"] * df["effGlo"] * df["effSta"]

        df["CHLT"] = 4 * (df["HLT2"]+df["HLT1"]+df["HLT0"]) * df["HLT2"] / (df["HLT1"] + 2*df["HLT2"])**2
                
        # ----- closure test for Z -> mumu
        # calculate true numbers at each step

        # from HLTTruth to RecoTruth
        df["NZ_reco"] = (df["HLT2"]+df["HLT1"]/2.)**2 / df["HLT2"]
        
        # from HLT to Truth
        df["NZ"] = df["NZ_reco"] / (df["effMu"])**2

        df["RZ"] = df["NZ"] / (df["Truth"])


        # ----- make plot for difference between reco and gen
        df["ReffID"]  = (df["effIDGen"]**2)  / (df["effID"]**2)
        df["ReffGlo"] = (df["effGloGen"]**2) / (df["effGlo"]**2)
        df["ReffSta"] = (df["effStaGen"]**2) / (df["effSta"]**2)
        df["RCHLT"]   = df["CHLT"] / df["CHLTGen"]
        
        df["RNZ_reco"] = df["NZ_reco"] / df["NZ_recoGen"] 

        ax = preparePlot("Reco / Gen", xVar)

        ax.plot([min(xx), max(xx+widthX)], [1,1], "k--")
        ax.plot(xx, df["RNZ_reco"].apply(lambda x: x.n).values, label="$N^\\mathrm{Z}_\\mathrm{reco}$", drawstyle='steps', color=colors[0])    
        ax.plot(xx, df["ReffSta"].apply(lambda x: x.n).values, label="$\\epsilon^{-2}_\\mathrm{Sta}$", drawstyle='steps', color=colors[2])    
        ax.plot(xx, df["ReffGlo"].apply(lambda x: x.n).values, label="$\\epsilon^{-2}_\\mathrm{Glo}$", drawstyle='steps', color=colors[1])    
        ax.plot(xx, df["ReffID"].apply(lambda x: x.n).values, label="$\\epsilon^{-2}_\\mathrm{ID}$", drawstyle='steps', color=colors[3])    
        # ax.plot(xx, df["RCHLT"].apply(lambda x: x.n).values, label="$C_\\mathrm{HLT}$", drawstyle='steps', color=colors[3])    

        finalizePlot(ax, "lower left", "GenReco", region="I", cmsTag="center left", year=year, nameTag=xVar, yRange=[0.981,1.024], ncol=2)   

        # ----- di muon acceptance
        print("# ----- Acceptance factor Gen:")

        acceptanceFactor = (df["NZGen"]).sum() / df["Truth"].sum()
        print("Acceptance factor (w/o correlation factors) {0}".format(acceptanceFactor))

        acceptanceFactor = (df["NZGen"] * df["CHLTGen"]).sum() / df["Truth"].sum()
        print("Acceptance factor for (w/ C_HLT ) {0}".format(acceptanceFactor))

        if useCorrIO:
            acceptanceFactor = (df["NZGen"] * df["CHLTGen"] * df["CID"]).sum() / df["Truth"].sum()
            print("Acceptance factor (w/ C_HLT C_ID) {0}".format(acceptanceFactor))

            # df["CIO"] = 1./df["CIO"]**0.5

            acceptanceFactor = (df["NZGen"] * df["CHLTGen"] * df["CID"] * df["CIO"]**2).sum() / df["Truth"].sum()
            print("Acceptance factor (w/ c_io^2 c_ID c_HLT) {0}".format(acceptanceFactor))


        print("# ----- Acceptance factor Reco:")

        acceptanceFactor = (df["NZ"]).sum() / df["Truth"].sum()
        print("Acceptance factor (w/o correlation factors) {0}".format(acceptanceFactor))

        acceptanceFactor = (df["NZ"] * df["CHLT"]).sum() / df["Truth"].sum()
        print("Acceptance factor for (w/ C_HLT ) {0}".format(acceptanceFactor))

        if useCorrIO:
            acceptanceFactor = (df["NZ"] * df["CHLT"] * df["CID"]).sum() / df["Truth"].sum()
            print("Acceptance factor (w/ C_HLT C_ID) {0}".format(acceptanceFactor))

            # df["CIO"] = 1./df["CIO"]**0.5

            acceptanceFactor = (df["NZ"] * df["CHLT"] * df["CID"] * df["CIO"]**2).sum() / df["Truth"].sum()
            print("Acceptance factor (w/ c_io^2 c_ID c_HLT) {0}".format(acceptanceFactor))



        if xVar == "nPV":
            # events with 0 PV have very low statistics and are set to the value of the first bin
            df.at[0,"RZ"] = df["RZ"].values[1]
            df.at[0,"CHLT"] = 1.0
        
        # if year[1] == "Summer16preVFP":
        #     rangeY=[0.82,1.025]
        # elif year[1] == "Summer16postVFP":
        #     rangeY=[0.95,1.025]
        # elif year[1] == "Fall17":
        #     rangeY = [0.98, 1.01]
        # else:
        # rangeY=[0.9875,1.0125]
        
        ax = preparePlot("Measurement / Truth", xVar)
        ax.plot([min(xx), max(xx+widthX)], [1,1], "k--")

        # yy0 = df["RZ"].apply(lambda x: x.n).values

        # yyerr0 = df["RZ"].apply(lambda x: x.s).values
        # xxerr0 = np.ones(len(xx))*widthX/2.

        # ax.errorbar(xx, yy0, yerr=yyerr0, xerr=xxerr0,
        #     # label="$N_\\mathrm{Z}\ C_\\mathrm{IO}^2 C_\\mathrm{ID} C_\\mathrm{HLT}$", 
        #     label="$N_\\mathrm{Z}$",
        #     #drawstyle='steps-mid', 
        #     color="black", 
        #     # color=colors[0], 
        #     linestyle="",marker=".", zorder=0)    

        yy0 = df["RZ"].apply(lambda x: x.n).values
        yy0 = np.append(yy0, yy0[-1])
        # ax.plot(xx0, yy0, label="$N^\\mathrm{Z}$", drawstyle='steps-mid', color=colors[1], zorder=0)                  

        # do_fits = False
        # if do_fits:
        #     # perform fit to quantify linearity
        #     func = linear
        #     popt, pcov = curve_fit(func, xx[1:], df["RZ_corr1"].apply(lambda x: x.n).values[1:])

        #     perr = np.sqrt(np.diag(pcov))
        #     # correlated values
        #     params = unc.correlated_values(popt, pcov)

        #     yy_corr1 = np.array([func(x, *params).n for x in xx])

        #     func = linear
        #     popt, pcov = curve_fit(func, xx[1:], df["RZ"].apply(lambda x: x.n).values[1:])

        #     perr = np.sqrt(np.diag(pcov))
        #     # correlated values
        #     params = unc.correlated_values(popt, pcov)

        #     yy = np.array([func(x, *params).n for x in xx])

        #     ax.plot(xx, yy, color=colors[1])    
        #     ax.plot(xx, yy_corr1, color=colors[2])      

        if useCorrIO:
            
            df["RZ_corrIO"] = df["RZ"] * df["CIO"]**2
            df["RZ_corrIOandID"] = df["RZ_corrIO"] * df["CID"]
            df["RZ_corrIOandIDandHLT"] = df["RZ_corrIOandID"] * df["CHLT"]

            yy0 = df["RZ_corrIO"].apply(lambda x: x.n).values
            yy0 = np.append(yy0, yy0[-1])

            # ax.plot(xx0, yy0, label="$N^\\mathrm{Z}\ C_\\mathrm{IO}^2$", drawstyle='steps-mid', color=colors[2], zorder=1)    
                
            yy0 = df["RZ_corrIOandID"].apply(lambda x: x.n).values
            yy0 = np.append(yy0, yy0[-1])

            # ax.plot(xx0, yy0, label="$N^\\mathrm{Z}\ C_\\mathrm{IO}^2 C_\\mathrm{ID}$", drawstyle='steps-mid', color=colors[3], zorder=2)            

            # xx0 = xx+widthX/2.
            yy0 = df["RZ_corrIOandIDandHLT"].apply(lambda x: x.n).values
            yy0 = np.append(yy0, yy0[-1])

            ax.plot(xx0, yy0, label="$N^\\mathrm{Z}\ C_\\mathrm{IO}^2 C_\\mathrm{ID} C_\\mathrm{HLT}$", drawstyle='steps-mid', color=colors[0], zorder=3)   

            # yyerr0 = df["RZ_corrIOandIDandHLT"].apply(lambda x: x.s).values
            # xxerr0 = np.ones(len(xx))*widthX/2.

            # ax.errorbar(xx0, yy0, yerr=yyerr0, xerr=xxerr0,
            #     label="$N_\\mathrm{Z}\ C_\\mathrm{IO}^2 C_\\mathrm{ID} C_\\mathrm{HLT}$", 
            #     # label="$N_\\mathrm{Z}$ corrected",
            #     #drawstyle='steps-mid', 
            #     #color="black", 
            #     color=colors[3], 
            #     linestyle="",marker=".", zorder=2)            

            # --> Do fit
            # print(">>> make fit")
            # func = linear #pol2 #quad
            # popt, pcov = curve_fit(func, xx0, yy0, sigma=yyerr0, absolute_sigma=True)

            # perr = np.sqrt(np.diag(pcov))
            # params = unc.correlated_values(popt, pcov)

            # f = lambda x: func(x, *params)
            # xFit = np.arange(0,100,0.5)
            # yFit = np.array([f(x).n for x in xFit])
            # yErrFit = np.array([f(x).s for x in xFit])

            # ax.plot(xFit, yFit, color=colors[0])   
            # <--

        else:
            df["RZ_corrHLT"] = df["RZ"] * df["CHLT"]

            # ax.errorbar(xx+widthX/2., df["RZ_corrHLT"].apply(lambda x: x.n).values, 
            #     yerr=df["RZ_corrHLT"].apply(lambda x: x.s).values, xerr=np.ones(len(xx))*widthX/2.,
            #     label="$N_\\mathrm{Z}\ C_\\mathrm{HLT}$", 
            #     #drawstyle='steps-mid', 
            #     color="black",#colors[0], 
            #     linestyle="",marker=".", zorder=2)              

        finalizePlot(ax, "lower left", "Reco", region="I", cmsTag="center left", year=year, nameTag=xVar, yRange=[0.991,1.029])                
            

        # save nonclosure as final correction for "kinematic selection"

        df["CAcceptance"] = 1./df["RZ_corrIOandIDandHLT"].apply(lambda x: x.n)

        outputname = "{0}/c_{1}_{2}.root".format(output, xVar, year[0])
        tFile = ROOT.TFile(outputname, "RECREATE")    

        for key in ["CID", "CIO", "CAcceptance"]:

            thc = ROOT.TH1D("{0}_{1}".format(key, "I"), 
                "{0} {1}".format(key,"I"), nBinsX, minX, maxX)
            for x, y in df[["x", key]].values:
                thc.SetBinContent(thc.FindBin(x), y)

            thc.Write()
            thc.Reset()
            thc.Delete()

