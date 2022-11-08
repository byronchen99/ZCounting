# make histograms of the parameters that are obtained in the fits

import ROOT
import glob
import pdb
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse
import os,sys
import json

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino",],
    "font.size": 16,
    'text.latex.preamble': [r"""\usepackage{bm}"""]
})

mpl.rcParams.update({
    "legend.fontsize" : "medium",
    "axes.labelsize" : "medium",
    "axes.titlesize" : "medium",
    "xtick.labelsize" : "medium",
    "ytick.labelsize" : "medium",
})

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", required=True, type=str, help="Give input dir")
parser.add_argument("-o","--output", default='./', type=str, help="Give output dir")
args = parser.parse_args()

outDir = args.output
if not os.path.isdir(outDir):
    os.mkdir(outDir)

inputDir = args.input

etaRegion = "I"
etaRegionZ = "I"

sig_mean_1_HLT = []
sig_sigma_1_HLT = []
sig_mean_2_HLT = []
sig_sigma_2_HLT = []

# alias: variable name
vars = {
    "sig_mean_HLT_1":       "sig_meanPass_1",
    "sig_sigma_HLT_1":      "sig_sigmaPass_1",
    "sig_mean_HLT_2":       "sig_meanPass_2",
    "sig_sigma_HLT_2":      "sig_sigmaPass_2",
    # "bkg_alpha_HLT_1":      "bkg_alphaPass_1",
    # "bkg_beta_HLT_1":       "bkg_betaPass_1",
    # "bkg_gamma_HLT_1":      "bkg_gammaPass_1",
    # "bkg_alpha_HLT_2":      "bkg_alphaPass_2",
    # "bkg_beta_HLT_2":       "bkg_betaPass_2",
    # "bkg_gamma_HLT_2":      "bkg_gammaPass_2",
    "bkg_t1_HLT_1":      "bkg_t1Pass_1",
    "bkg_t1_HLT_2":      "bkg_t1Pass_2",
    # "sig_meanPass_Sel":     "sig_meanPass_0",
    # "sig_sigmaPass_Sel":    "sig_sigmaPass_0",
    "sig_meanFail_Sel":     "sig_meanFail_0",
    "sig_sigmaFail_Sel":    "sig_sigmaFail_0",
    # "bkg_alphaPass_Sel":      "bkg_alphaPass_0",
    # "bkg_betaPass_Sel":       "bkg_betaPass_0",
    # "bkg_gammaPass_Sel":      "bkg_gammaPass_0",
    "bkg_alphaFail_Sel":      "bkg_alphaFail_0",
    "bkg_betaFail_Sel":       "bkg_betaFail_0",
    "bkg_gammaFail_Sel":      "bkg_gammaFail_0",    
    # "sig_meanPass_Glo":     "sig_meanPass_0",
    # "sig_sigmaPass_Glo":    "sig_sigmaPass_0",
    # "sig_meanFail_Glo":     "sig_meanFail_0",
    # "sig_sigmaFail_Glo":    "sig_sigmaFail_0",
    # "bkg_alphaPass_Glo":      "bkg_alphaPass_0",
    # "bkg_betaPass_Glo":       "bkg_betaPass_0",
    # "bkg_gammaPass_Glo":      "bkg_gammaPass_0",
    # "bkg_alphaFail_Glo":      "bkg_alphaFail_0",
    # "bkg_betaFail_Glo":       "bkg_betaFail_0",
    # "bkg_gammaFail_Glo":      "bkg_gammaFail_0",    
    # "sig_meanPass_Sta":     "sig_meanPass_0",
    # "sig_sigmaPass_Sta":    "sig_sigmaPass_0",
    # "sig_meanFail_Sta":     "sig_meanFail_0",
    # "sig_sigmaFail_Sta" :   "sig_sigmaFail_0",
    # "bkg_alphaPass_Sta":      "bkg_alphaPass_0",
    # "bkg_betaPass_Sta":       "bkg_betaPass_0",
    # "bkg_gammaPass_Sta":      "bkg_gammaPass_0",
    # "bkg_alphaFail_Sta":      "bkg_alphaFail_0",
    # "bkg_betaFail_Sta":       "bkg_betaFail_0",
    # "bkg_gammaFail_Sta":      "bkg_gammaFail_0",
    "sig_meanPass_Trk":     "sig_meanPass_0",
    "sig_sigmaPass_Trk":    "sig_sigmaPass_0",
    "sig_meanFail_Trk":     "sig_meanFail_0",
    "sig_sigmaFail_Trk" :   "sig_sigmaFail_0",
    # "bkg_alphaPass_Trk":      "bkg_alphaPass_0",
    # "bkg_betaPass_Trk":       "bkg_betaPass_0",
    # "bkg_gammaPass_Trk":      "bkg_gammaPass_0",
    "bkg_t1Pass_Trk":      "bkg_t1Pass_0",
    "bkg_alphaFail_Trk":      "bkg_alphaFail_0",
    "bkg_betaFail_Trk":       "bkg_betaFail_0",
    "bkg_gammaFail_Trk":      "bkg_gammaFail_0",
    "sig_mean_Trk_2":       "sig_meanPass_0",
    "sig_sigma_Trk_2":      "sig_sigmaPass_0",
    "bkg_alpha_Trk_2":      "bkg_alphaPass_0",
    "bkg_beta_Trk_2":       "bkg_betaPass_0",
    "bkg_gamma_Trk_2":      "bkg_gammaPass_0",
}

xlabels = {
    "mean" : "mean",
    "sigma": "sigma",
    "alpha" : "alpha",
    "beta" : "beta",
    "gamma" : "gamma",
    "t1" : "t",
    "statusMinuit": "Minuit status",
    "statusHesse": "Hesse status",
    "statusMinos": "Minos status",
    "covQual": "Covariance matrix quality"
}

xlimits = {
    "mean" : [-2.5,2.5],
    "sigma": [0.1, 5],
    "alpha" : [0, 300],
    "beta" : [0.0, 0.1],
    "gamma" : [0.0, 0.1],
    "t1" : [-1.0, 0.0],
    "statusMinuit": [0,5],
    "statusHesse": [0,4],
    "statusMinos": [0,4],
    "covQual": [0,3]
}


def plot_variable(key, xx, era="RunII"):

    mean = xx.mean()
    std = xx.std()

    plt.clf()
    # fig = plt.figure(figsize=(6.4,8.))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    fig.subplots_adjust(hspace=0.0, left=0.13, right=0.98, top=0.99, bottom=0.13)
    
    if len(key.split("_")) > 2:
        var_name = key.split("_")[1].split("Fail")[0].split("Pass")[0]

        labelname = " ".join(key.split("_")[2:])
        if "Fail" in key:
            labelname += " Fail"
        elif "Pass" in key:
            labelname += " Pass"
    else:
        var_name = key.split("_")[0]
        labelname = key.split("_")[1]

    ax.set_ylabel("Frequency")
    ax.set_xlabel(xlabels[var_name])

    nEntries, bins, _ = ax.hist(xx, bins=100, range=xlimits[var_name], label=labelname)

    ax.text(0.03, 0.97, "\\bf{CMS}", verticalalignment='top', transform=ax.transAxes, weight="bold")
    ax.text(0.14, 0.97, "\\emph{Work in progress}", verticalalignment='top', transform=ax.transAxes,style='italic')
    ax.text(0.03, 0.91, era, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes,style='italic')
    ax.text(0.98, 0.85, "$\\mu$ = {0}".format(round(mean,3)), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,style='italic')
    ax.text(0.98, 0.78, "$\\sigma$ = {0}".format(round(std,3)), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,style='italic')


    leg = plt.legend(loc="upper right")

    # leg = ax1.legend([(p2[0], p1[0]), p3], ['Simulation', 'Data'], loc="upper left", frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
    # 
    leg.get_frame().set_linewidth(0.8)

    # ax1.yaxis.set_label_coords(-0.12, 0.5)

    ax.set_xlim(xlimits[var_name])
    ax.set_ylim([0, max(nEntries)*1.2])

    plt.savefig("{0}/fits_{1}_{2}.png".format(outDir, key, era))
    plt.savefig("{0}/fits_{1}_{2}.pdf".format(outDir, key, era))
    plt.close()

if args.input.endswith("json"):
    with open(args.input, "r") as file_info:
        info = json.load(file_info)

    for key, var in info.items():
        
        xx = np.array(var)
        if len(xx) == 0:
            continue
        
        plot_variable(key, xx)


    exit()

for iEra, run_lo, run_hi in (
    ("2016preVFP", 272007, 278808),
    ("2016postVFP", 278809, 294645),
    ("2017", 297020, 307083),
    ("2018", 315252, 325273),
    # ("2022", 355100, 999999),
):
    print("Now at era {0}".format(iEra))

    # set empty lists
    info = {}
    for var in vars:
        info[var] = []

    for var in ("covQual", "statusHesse", "statusMinos", "statusMinuit"):
        for eff in ("HLT", "Sel", "Trk"):
            info["{0}_{1}".format(var,eff)] = []
   

    for dir_run in glob.glob(inputDir+"/Run*"):
        run = int(dir_run.split("Run")[-1].split("to")[0])
        if run < run_lo or run > run_hi:
            continue
        print("Now at Run {0}".format(run))
        # from HLT efficiency and yield fits
        for file_m in glob.glob(dir_run+"/workspace_{0}_*.root".format(etaRegionZ)):
            f = ROOT.TFile(file_m,"READ")
            w = f.Get("workspace")
            
            for key, var in filter(lambda x: "HLT" in x[0], vars.items()):
                info[key].append(w.var(var).getVal())

            # check fit result for fit quality
            r = f.Get("fitResult")
            info["covQual_HLT"].append(r.covQual())
            info["statusHesse_HLT"].append(r.status()//100)
            info["statusMinos_HLT"].append((r.status()//10)%10)
            info["statusMinuit_HLT"].append(r.status()%10)

            f.Close()
            f.Delete()
            w.Delete()
            r.Delete()

        # from efficiency fits
        for eff in ["Trk",]:
            for file_m in glob.glob(dir_run+"/workspace_{0}_{1}_*.root".format(eff, etaRegion)):
                f = ROOT.TFile(file_m,"READ")
                w = f.Get("workspace")
                
                for key, var in filter(lambda x: x[0].endswith(eff), vars.items()):
                    info[key].append(w.var(var).getVal())

                # check fit result for fit quality
                r = f.Get("fitResult")
                info["covQual_{0}".format(eff)].append(r.covQual())
                info["statusHesse_{0}".format(eff)].append(r.status()//100)
                info["statusMinos_{0}".format(eff)].append((r.status()//10)%10)
                info["statusMinuit_{0}".format(eff)].append(r.status()%10)
                f.Close()
                f.Delete()
                w.Delete()
                r.Delete()

        # from yield fits
        for eff in ["Trk_{0}_2".format(etaRegion), "Sel_{0}".format(etaRegion),]: 
            for file_m in glob.glob(dir_run+"/workspace_yield_{0}_*.root".format(eff)):
                f = ROOT.TFile(file_m,"READ")
                w = f.Get("workspace")

                #string for comparison
                vareff = eff.replace("_{0}".format(etaRegion),"")
                
                for key, var in filter(lambda x: x[0].endswith(vareff), vars.items()):
                    info[key].append(w.var(var).getVal())

                if "Trk" in eff:
                    f.Close()
                    f.Delete()
                    w.Delete()
                    continue
                    
                # check fit result for fit quality
                r = f.Get("fitResult")
                info["covQual_Sel"].append(r.covQual())
                info["statusHesse_Sel"].append(r.status()//100)
                info["statusMinos_Sel"].append((r.status()//10)%10)
                info["statusMinuit_Sel"].append(r.status()%10)

                f.Close()
                f.Delete()
                w.Delete()
                r.Delete()


    with open(outDir+"/info_{0}.json".format(iEra),"w") as outfile:
        json.dump(info, outfile, indent=4, sort_keys=True)

    for key, var in info.items():
        xx = np.array(var)
        if len(xx) == 0:
            continue
        
        plot_variable(key, xx, era=iEra)

