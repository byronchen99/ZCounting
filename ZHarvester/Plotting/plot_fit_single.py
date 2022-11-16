import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

import json
import numpy as np
import uncertainties as unc
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import pdb
import os
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

textsize = 16
markersize = 4
xlabel = "Number of primary vertices"
ylabel = "Correlation factor"

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino",],
    "font.size": textsize,
    'text.latex.preamble': [r"""\usepackage{bm}"""],
    'figure.figsize': [7, 4.5]
})

mpl.rcParams.update({
    "legend.fontsize" : "medium",
    "axes.labelsize" : "medium",
    "axes.titlesize" : "medium",
    "xtick.labelsize" : "medium",
    "ytick.labelsize" : "medium",
})

xlabel = "Z candidate mass [GeV]"
markersize = 4

parser = argparse.ArgumentParser()

parser.add_argument("-w","--workspace", required=True, type=str, help="Workspace in .root format")
parser.add_argument("--label",  default='Work in progress',  type=str, help="give label ")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

f1 = ROOT.TFile(args.workspace,"READ")

if not os.path.isdir(args.saveDir):
    os.mkdir(args.saveDir)

logscale = False
pulls = True

parts = args.workspace.split("/")[-1].split("_")[2:]

etaRegion = parts[1]
passRegion = int(parts[2])
efficiency = parts[0]
measurement = parts[3][:-5]

name = f"{efficiency}_{passRegion}"

passString = "Pass" if passRegion > 0 else "Fail"

# load RooCMS shape for fitting
ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__))+"/RooCMSShape.cc")

fitResult = f1.Get("fitResult")
ws = f1.Get("workspace")

# load RooDataHist, get each data set and convert to root histograms
data = ws.data("ZReco")

hist = data.createHistogram("m")

# get functions
mBkg = ws.pdf("background{0}_0".format(passString))
mSig = ws.pdf("signal{0}_0".format(passString))
# mTot = ws.pdf("totalPdf")

# get parameters
mass = ws.arg("m")

nBkg = ws.arg("Nbkg")
nSig = ws.arg("Nsig")
nTotal = nSig.getVal() + nBkg.getVal()

# fill arrays with function values

nBinsData = hist.GetNbinsX()
binWidthData = hist.GetBinWidth(0)

xMin = hist.GetBinLowEdge(1)
xMax = hist.GetBinLowEdge(nBinsData)+binWidthData

fData = np.array([hist.GetBinContent(i) for i in range(1, nBinsData+1)]) / binWidthData
fDataErr = np.array([hist.GetBinError(i) for i in range(1, nBinsData+1)]) / binWidthData
xData = np.arange(xMin+binWidthData/2., xMax+binWidthData/2., binWidthData)

# fTot = []
fSig = []
fBkg = []
for x in xData:
    mass.setVal(x)
    # fTot.append(mTot.evaluate())
    fSig.append(mSig.getValV())
    fBkg.append(mBkg.getValV())
    
fBkg = np.array(fBkg) / sum(fBkg) * nBkg.getVal() / binWidthData
fSig = np.array(fSig) / sum(fSig) * nSig.getVal() / binWidthData
# fTot = np.array(fTot) / sum(fTot) * nTotal / binWidthData
fTot = fBkg + fSig

# --- plot
plt.clf()

if pulls:        
    fig = plt.figure(figsize=(6.4,5.))
    gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    fig.subplots_adjust(left=0.15, right=0.97, top=0.99, bottom=0.125, hspace=0.0)
else:
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    fig.subplots_adjust(hspace=0.0, left=0.12, right=0.97, top=0.99, bottom=0.125)

# index = 2 if name == "Pass" else 1

if binWidthData == 1:
    ylabel = "Z boson candidates / 1 GeV"
else:
    ylabel = "Z boson candidates / {0} GeV".format(binWidthData)

ax1.set_xlabel(xlabel)
ax1.set_ylabel(ylabel)

xpos = 0.03
if pulls:
    pos = [0.48, 0.60, 0.72, 0.86, 0.95]
else:
    pos = [0.48, 0.60, 0.72, 0.86, 0.95]

ax1.text(xpos, pos[2], "$N^\mathrm{sig}_"+str(int(passRegion))+"$", verticalalignment='top', transform=ax1.transAxes)
ax1.text(xpos, pos[1], "$N^\mathrm{bkg}_"+str(int(passRegion))+"$", verticalalignment='top', transform=ax1.transAxes)    
ax1.text(xpos+0.09, pos[2], "$ = "+str(round(nSig.getVal(),1))+"\pm"+str(round(nSig.getPropagatedError(fitResult),1))+"$", verticalalignment='top', transform=ax1.transAxes)
ax1.text(xpos+0.09, pos[1], "$ = "+str(round(nBkg.getVal(),1))+"\pm"+str(round(nBkg.getError(),1))+"$", verticalalignment='top', transform=ax1.transAxes)

# if name == "Fail":
#     ax1.text(xpos, pos[0], "$\\epsilon_\mathrm{HLT}^{\mu}$", verticalalignment='top', transform=ax1.transAxes)
#     ax1.text(xpos+0.09, pos[0], "$ = "+str(round(eff.getVal(),3))+"\pm"+str(round(eff.getError(),3))+"$",verticalalignment='top', transform=ax1.transAxes)

# if name == "Pass":
#     ax1.text(xpos, pos[0], "$N^\mathrm{Z}_\mathrm{reco}$", verticalalignment='top', transform=ax1.transAxes)
#     # ax1.text(0.41, pos[0], "$N^\mathrm{Z} \\epsilon_\mathrm{ID}^{\mathrm{Z}}$", verticalalignment='top', transform=ax1.transAxes)
#     ax1.text(xpos+0.09, pos[0], "$ = "+str(round(nSigTotal.getVal(),1))+"\pm"+str(round(nSigTotal.getError(),1))+"$",verticalalignment='top', transform=ax1.transAxes)
    

# ax1.text(0.44, 0.43, "$N^\mathrm{sig}_"+str(int(index))+" = "+str(round(nSig.getVal(),1))+"\pm"+str(round(nSig.getPropagatedError(fitResult),1))+"$", verticalalignment='top', transform=ax1.transAxes)
# ax1.text(0.44, 0.34, "$N^\mathrm{bkg}_"+str(int(index))+" = "+str(round(nBkg.getVal(),1))+"\pm"+str(round(nBkg.getError(),1))+"$", verticalalignment='top', transform=ax1.transAxes)
# ax1.text(0.44, 0.25, "$\\epsilon_\mathrm{HLT}^{\mu} = "+str(round(eff.getVal(),3))+"\pm"+str(round(eff.getError(),3))+"$",
#     # +"^{+"+str(round(eff.getErrorHi(),3))+"}_{"+str(round(eff.getErrorLo(),3))+"}$"
#     verticalalignment='top', transform=ax1.transAxes)
ax1.text(xpos, pos[3], "$20\,\mathrm{pb}^{-1}\ \mathrm{at}\ \sqrt{s} = 13\,\mathrm{TeV}$",verticalalignment='top', transform=ax1.transAxes)
ax1.text(xpos, pos[4], "\\bf{CMS}", verticalalignment='top', transform=ax1.transAxes, weight="bold")
ax1.text(xpos+0.11, pos[4], "\\emph{"+args.label+"}", verticalalignment='top', transform=ax1.transAxes,style='italic')

ax1.step(xData, fTot, label="Sig.+Bkg.", color="red", zorder=2, where="mid")
ax1.step(xData, fBkg, label="Bkg.", color="blue", zorder=1, where="mid")
# ax1.plot(np.arange(69.5,250.5,1), fData, label="Data", marker="o", linewidth=0, color="black", zorder=3, markersize=markersize)
ax1.errorbar(xData, fData, yerr=fDataErr, label="Data", zorder=3,
    fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize)

leg = ax1.legend(loc="upper right",ncol=1, markerscale=1, frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
leg.get_frame().set_linewidth(0.8)

yMax = max(max(fTot),max(fData+fDataErr))

if logscale:
    ax1.set_ylim([0.1, yMax*1.5])
    ax1.set_yscale('log')
else:
    ax1.set_ylim([-0.1, yMax*1.025])
    
ax1.set_xlim([xMin, xMax])
#ax1.minorticks_on()
#ax1.tick_params(axis='both', which='both', direction="in", right=True, top=True)

if pulls:
    ax1.xaxis.set_major_locator(ticker.NullLocator())
    # ax1.set(xticklabels=[])
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("Pulls")
    pullErr = np.array([err if err > 0 else np.sqrt(fTot[i]) for i, err in enumerate(fDataErr)])
    yPulls = (fData - fTot) / pullErr 
    # ax2.errorbar(xData, yPulls, yerr=np.ones(len(yPulls)), zorder=3,
    #     fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize)
    ax2.plot(np.array([xMin, xMax]), np.array([0., 0.]), color="black",linestyle="-", linewidth=1)

    ax2.plot(xData, yPulls, "ko", markersize=markersize)
    ax2.set_xlim([xMin, xMax])
    ax2.set_ylim([-3., 3.])
    ax2.set_yticks([-2,0,2])
    ax1.yaxis.set_label_coords(-0.12, 0.5)
    ax2.yaxis.set_label_coords(-0.12, 0.5)

#ax1.yaxis.set_label_coords(-0.12, 0.5)
print("Plot {0}/fit_{1}.pdf".format(args.saveDir, name))
plt.savefig(args.saveDir+"/fit_{0}.pdf".format(name))
plt.close()
