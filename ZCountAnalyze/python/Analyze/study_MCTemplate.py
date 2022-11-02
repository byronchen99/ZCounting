#!/usr/bin/env python3
#
################################################################################
#
# study MC templates:
#   - impact of pileup (on mean, width, ...)
#
################################################################################

import uproot 
import awkward as ak
import numpy as np
import pdb
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import uncertainties as unc
from scipy.optimize import curve_fit

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import linear

import argparse
parser = argparse.ArgumentParser(prog='./histogramize')
parser.add_argument(
    '-i', '--input', nargs='+',
    help='specify input ntuple root files'
)
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output directory'
)
args = parser.parse_args()

textsize = 16
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

### acceptance cuts
ptCut = 25
etaCut = 2.4

mass_lo = 81
mass_hi = 101
mass_nBins = int(mass_hi-mass_lo)

pu_lo = 0
pu_hi = 75
pu_nBins = int(pu_hi-pu_lo)

pileup = "nPU"

xlabels = {
    "nPU":"Number of pileup", 
    "nPV":"Number of primary vertices"}
ylabels = {
    "avg":"mean(mass)",
    "std":"std(mass)",
    "mass": "mass"}

### resources

dirOut = args.output

if not os.path.isdir(dirOut):
    print(f"create output directory {dirOut}")
    os.mkdir(dirOut)
    

# branches = ["mass", "ptTag", "ptProbe", "etaTag", "etaProbe", "nPV","pass","match1","match2","eventWeight"]
branches = ["mass", pileup, "pass"]
aliases = {

}


hists = dict()
for treename in ["HLT", "Sel", "Trk"]:
    print(f"Now at {treename}")    

    filenames = [f+":"+treename for f in args.input]

    # process the events in batches to avoid using too much memory
    for batch, report in uproot.iterate(filenames, step_size="100 MB", library="ak", 
        filter_name=branches, aliases=aliases, report=True,
        cut=f"(mass > {mass_lo}) & (mass < {mass_hi}) \
            & (ptTag > {ptCut}) & (ptProbe > {ptCut}) \
            & (abs(etaTag) < {etaCut}) & (abs(etaProbe) < {etaCut}) \
            & match1 & match2"
    ):
        if len(batch) == 0:
            continue
        print(f"Have a batch with {len(batch)} events")    

        year = report.file_path.split("-")[-2]

        if treename=="HLT":
            filters = {
                "0": batch["pass"] == 0,
                "1": batch["pass"] == 1,
                "2": batch["pass"] == 2
                }
        elif treename=="Sel":
            filters = {
                "fail": batch["pass"] == 0}      
        else:
            filters = {
                "pass": batch["pass"] >= 1,
                "fail": batch["pass"] == 0}
        
        # loop over pair collections with names to store them in dictionary
        for name, filter in filters.items():

            xx = batch[filter]["mass"].to_numpy()
            yy = batch[filter][pileup].to_numpy()

            avgs = np.array(
                [np.median(batch[filter & (batch[pileup] == i)]["mass"]) for i in range(pu_lo, pu_hi)])


            stds = np.array(
                [np.std(batch[filter & (batch[pileup] == i)]["mass"]) for i in range(pu_lo, pu_hi)])

            histo = np.histogram2d(xx, yy, bins=[mass_nBins, pu_nBins], 
                range=((mass_lo, mass_hi), (pu_lo, pu_hi)))

            h_std = stds, histo[2]
            h_avg = avgs, histo[2]

            histname = f"h2d_mass_{pileup}_{treename}_{name}_{year}"
            hstdname = f"h1d_std_{pileup}_{treename}_{name}_{year}"
            havgname = f"h1d_avg_{pileup}_{treename}_{name}_{year}"

            if histname in hists.keys():
                if len(hists[histname]) == 3:
                    hists[histname] = (histo[0] + hists[histname][0], histo[1], histo[2])
                else:
                    hists[histname] = (histo[0] + histo[histname][0], histo[1])
            else:
                hists[histname] = histo
                hists[hstdname] = h_std
                hists[havgname] = h_avg

print("Make plots")
# plot histograms
for name, histo in hists.items():
    if len(histo) == 3:
        #2d histogram
        continue
    elif len(histo) == 2:
        #1d histogram
        plt.clf()
        fig = plt.figure()
        fig.subplots_adjust(left=0.15, right=0.98, top=0.97, bottom=0.125)
        ax = fig.add_subplot(111)

        xT = 0.02
        yT = 0.98

        xx = (histo[1][1:] - histo[1][:-1])/2. + histo[1][:-1]
        yy = histo[0]

        ax.plot(xx, yy, "k.", label=name)

        func = linear #pol2 #quad
        popt, pcov = curve_fit(func, xx, yy)

        perr = np.sqrt(np.diag(pcov))
        params = unc.correlated_values(popt, pcov)

        f = lambda x: func(x, *params)
        yFit = np.array([f(x).n for x in xx])
        yFitErr = np.array([f(x).s for x in xx])

        ax.plot(xx, yFit, "k--", label=name)

        ax.fill_between(xx, yFit - yFitErr, yFit + yFitErr,
                    color='grey', alpha=0.2, zorder=1) 

        ax.set_xlabel(xlabels[name.split("_")[2]])
        ax.set_ylabel(ylabels[name.split("_")[1]])

        ax.text(xT, yT, "\\bf{CMS}", verticalalignment='top', transform=ax.transAxes, weight="bold")
        ax.text(xT+0.11, yT, "\\emph{Simulation}", verticalalignment='top', transform=ax.transAxes,style='italic')
        ax.text(xT, yT-0.07, "\\emph{Work in progress}", verticalalignment='top', transform=ax.transAxes,style='italic')    
        ax.text(0.98, yT, yearnames[name.split("_")[-1]], verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)
        ax.text(0.98, yT-0.07, " ".join(name.split("_")[3:-1]) , verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)

        plt.savefig(f"{dirOut}/{name}.png")
        plt.savefig(f"{dirOut}/{name}.pdf")



# open output root file
print(f"Write out results in `output.root`")
output = uproot.recreate(f"{dirOut}/output.root")

# write the histograms out
for name, histo in hists.items():
    output[name] = histo




