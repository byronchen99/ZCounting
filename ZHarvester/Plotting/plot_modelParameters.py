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

from common import parsing, plotting, logging

parser = parsing.parser_plot()
parser.add_argument("-i","--input", required=True, type=str, help="Give input dir")
args = parser.parse_args()
log = logging.setup_logger(__file__, args.verbose)

outDir = args.output
if not os.path.isdir(outDir):
    os.mkdir(outDir)

inputDir = args.input

colors, textsize, labelsize, markersize = plotting.set_matplotlib_style()


etaRegion = "I"
etaRegionZ = "I"


xlabels = {
    "mean" : "mean",
    "sigma": "sigma",
    "n" : "n",
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
    "n" : [0.5, 10],
    "alpha" : [0, 80],
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
        var_name = key.split("_")[-1].split("Fail")[0].split("Pass")[0]

        labelname = key.split("_")[0]
        if "Fail" in key:
            labelname += " Fail"
        elif "HLT_1" in key:
            labelname += " 1"
        elif "HLT_2" in key:
            labelname += " 2"
        elif "Pass" in key:
            labelname += " Pass"

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

    iEra = args.input.split(".")[-2].split("_")[-1]

    for key, var in info.items():
        
        xx = np.array(var)
        if len(xx) == 0:
            continue
        
        plot_variable(key, xx, era=iEra)


    exit()

for iEra, run_lo, run_hi in (
    ("2016preVFP", 272007, 278808),
    ("2016postVFP", 278809, 294645),
    ("2017", 297020, 306925),
#    ("2017H", 306926, 307083),
    ("2018", 315252, 325273),
    # ("2022", 355100, 999999),
):
    print("Now at era {0}".format(iEra))

    # set empty lists
    info = {}

    for var in ("covQual", "statusHesse", "statusMinos", "statusMinuit"):
        for eff in ("HLT", "ID", "Glo", "Sta"):
            for ipass in [0, 1, 2]:

                if eff == "HLT":
                    pass_key = str(ipass)
                elif ipass == 0:
                    pass_key = "Fail"
                elif ipass == 1:
                    pass_key = "Pass"
                else: 
                    continue      

                info[f"{eff}_{pass_key}_{var}"] = []
   

    for dir_run in glob.glob(inputDir+"/Run*"):
        run = int(dir_run.split("Run")[-1].split("to")[0])
        if run < run_lo or run > run_hi:
            continue
        if run == 275310:
            continue
        print("Now at Run {0}".format(run))

        # # from HLT efficiency and yield fits
        # for file_m in glob.glob(dir_run+"/workspace_{0}_*.root".format(etaRegionZ)):
        #     f = ROOT.TFile(file_m,"READ")
        #     w = f.Get("workspace")
            
        #     for key, var in filter(lambda x: "HLT" in x[0], vars.items()):
        #         info[key].append(w.var(var).getVal())

        #     # check fit result for fit quality
        #     r = f.Get("fitResult")
        #     info["covQual_HLT"].append(r.covQual())
        #     info["statusHesse_HLT"].append(r.status()//100)
        #     info["statusMinos_HLT"].append((r.status()//10)%10)
        #     info["statusMinuit_HLT"].append(r.status()%10)

        #     f.Close()
        #     f.Delete()
        #     w.Delete()
        #     r.Delete()

        # # from efficiency fits
        # for eff in ["Trk",]:
        #     for file_m in glob.glob(dir_run+"/workspace_{0}_{1}_*.root".format(eff, etaRegion)):
        #         f = ROOT.TFile(file_m,"READ")
        #         w = f.Get("workspace")
                
        #         for key, var in filter(lambda x: x[0].endswith(eff), vars.items()):
        #             info[key].append(w.var(var).getVal())

        #         # check fit result for fit quality
        #         r = f.Get("fitResult")
        #         info["covQual_{0}".format(eff)].append(r.covQual())
        #         info["statusHesse_{0}".format(eff)].append(r.status()//100)
        #         info["statusMinos_{0}".format(eff)].append((r.status()//10)%10)
        #         info["statusMinuit_{0}".format(eff)].append(r.status()%10)
        #         f.Close()
        #         f.Delete()
        #         w.Delete()
        #         r.Delete()

        # from fits of single histograms
        for eff in ["HLT","ID","Glo","Sta"]: 
            for ipass in [0, 1, 2]:

                if eff == "HLT":
                    pass_key = str(ipass)
                elif ipass == 0:
                    pass_key = "Fail"
                elif ipass == 1:
                    pass_key = "Pass"
                else: 
                    continue                


                for file_m in glob.glob(dir_run+"/workspace_yield_{0}_{1}_{2}*.root".format(eff,etaRegion,ipass)):
                    f = ROOT.TFile(file_m,"READ")
                    w = f.Get("workspace")

                    for p in w.allVars():

                        title = p.GetName().replace("Fail","").replace("_0","")
                        if "chi" in title or "peak" in title:
                            continue
                        if title in ["Nsig", "Nbkg", "m", "bkg_peak"]:
                            continue

                        key = f"{eff}_{pass_key}_{title}"
                        if key not in info.keys():
                            info[key] = [p.getVal(),]
                        else:
                            info[key].append(p.getVal())

                    # check fit result for fit quality
                    r = f.Get("fitResult")
                    info[f"{eff}_{pass_key}_covQual"].append(r.covQual())
                    info[f"{eff}_{pass_key}_statusHesse"].append(r.status()//100)
                    info[f"{eff}_{pass_key}_statusMinos"].append((r.status()//10)%10)
                    info[f"{eff}_{pass_key}_statusMinuit"].append(r.status()%10)

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

