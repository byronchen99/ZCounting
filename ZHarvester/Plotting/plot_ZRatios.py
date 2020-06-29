import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas as pd
import numpy as np

import pdb

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("--rates", required=True, type=str, help="csv file with z rates")
parser.add_argument("--xsec",  required=True, type=str, help="csv file where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)


def make_ratio(data_lo, data_hi, run_range=None, name=""):

    if run_range:
        data_hi = data_hi.query("run >= {0} & run < {1}".format(*run_range))

    rLumi = data_hi['lumiRec'].sum() / data_lo['lumiRec'].sum()

    rZ_BB = data_hi['zDelBB_mc'].sum() / data_lo['zDelBB_mc'].sum()
    rZ_BE = data_hi['zDelBE_mc'].sum() / data_lo['zDelBE_mc'].sum()
    rZ_EE = data_hi['zDelEE_mc'].sum() / data_lo['zDelEE_mc'].sum()
    rZ_tot = (data_hi['zDelBB_mc'].sum() + data_hi['zDelBE_mc'].sum() + data_hi['zDelEE_mc'].sum()) / (data_lo['zDelBB_mc'].sum() + data_lo['zDelBE_mc'].sum() + data_lo['zDelEE_mc'].sum())

    rZ_BB_err = rZ_BB * 1./np.sqrt(data_xsec['zYieldBB'].sum())
    rZ_BE_err = rZ_BE * 1./np.sqrt(data_xsec['zYieldBE'].sum())
    rZ_EE_err = rZ_EE * 1./np.sqrt(data_xsec['zYieldEE'].sum())
    rZ_tot_err = rZ_tot / (data_lo['zDelBB_mc'].sum() + data_lo['zDelBE_mc'].sum() + data_lo['zDelEE_mc'].sum()) * np.sqrt(data_lo['zYieldBB'].sum() + data_lo['zYieldBE'].sum() + data_lo['zYieldEE'].sum())

    points = np.array([rLumi, rZ_BB, rZ_BE, rZ_EE, rZ_tot])
    points_err = np.array([rLumi*0.013, rZ_BB_err, rZ_BE_err, rZ_EE_err, rZ_tot_err])


    ########## Plot ##########
    xmin = 0.9
    xmax = 1.5

    graphs = []
    for i, (ipoint, ierr, nam, ptr) in enumerate(
        zip(points,points_err,
            ['Lumi','Z BB','Z BE','Z EE','Z total'],
            [20, 21, 22, 23, 34],
        )):
        graph=ROOT.TGraphErrors(1, np.array([0.1*i+1]), np.array([ipoint]),  np.array([0.]), np.array(ierr))
        graph.SetName(nam)
        graph.SetTitle(nam)
        graph.SetMarkerStyle(ptr)
        graph.SetMarkerColor(i+1)
        graph.SetFillStyle(1001)
        graph.SetMarkerSize(1.5)
        graphs.append(graph)

    c2=ROOT.TCanvas("c2","c2",500,600)
    pad1 = ROOT.TPad("pad1", "pad1", 0., 0.4, 1, 1.0)
    pad1.SetBottomMargin(0.)
    c2.SetTicks()
    pad1.SetLeftMargin(0.2)
    pad1.SetRightMargin(0.01)
    pad1.SetTopMargin(0.1)
    pad1.SetTickx()
    pad1.SetTicky()
    pad1.Draw()
    pad1.cd()

    textsize = 24./(pad1.GetWh()*pad1.GetAbsHNDC())

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(11)
    latex.SetTextFont(42)
    latex.SetTextSize(textsize)

    ymin = min(points-points_err)*0.999
    ymax = ymin + 1.15 * (max(points) - ymin)

    graphs[0].GetYaxis().SetRangeUser(ymin, ymax)
    graphs[0].GetXaxis().SetRangeUser(xmin, xmax)
    graphs[0].GetXaxis().SetLabelSize(0)
    graphs[0].GetYaxis().SetTitle("Ratio")
    graphs[0].GetYaxis().SetTitleOffset(1.2)
    graphs[0].GetYaxis().SetTitleSize(textsize)
    graphs[0].GetYaxis().SetLabelSize(textsize)
    graphs[0].Draw("AP")

    legend=ROOT.TLegend(0.75,0.55,0.98,0.85)

    latex.DrawLatex(0.2, 0.91, name)


    for graph in graphs:
        graph.Draw("P same")
        legend.AddEntry(graph,"","pe")

    legend.SetTextFont(42)
    legend.SetTextSize(textsize)

    legend.Draw("same")

    graphs[0].SetTitle("")
    graphs[0].Draw("same")

    ### ratio ###

    points_err /= points[0]
    points /= points[0]

    rgraphs = []
    for i, (ipoint, ierr, nam, ptr) in enumerate(
        zip(points,points_err,
            ['Lumi','Z BB','Z BE','Z EE','Z total'],
            [20, 21, 22, 23, 34],
        )):
        graph=ROOT.TGraphErrors(1, np.array([0.1*i+1]), np.array([ipoint]),  np.array([0.]), np.array(ierr))
        graph.SetName(nam)
        graph.SetTitle(nam)
        graph.SetMarkerStyle(ptr)
        graph.SetMarkerColor(i+1)
        graph.SetFillStyle(1001)
        graph.SetMarkerSize(1.5)
        rgraphs.append(graph)

    c2.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.4)
    pad2.SetLeftMargin(0.2)
    pad2.SetRightMargin(0.01)
    pad2.SetTopMargin(0.0)
    pad2.SetBottomMargin(0.001)
    pad2.SetTickx()
    pad2.SetTicky()
    pad2.Draw("ALPF")
    pad2.cd()

    textsize = 24./(pad2.GetWh()*pad2.GetAbsHNDC())

    ymin = min(points-points_err)*0.999
    ymax = ymin + 1.15 * (max(points) - ymin)

    rgraphs[0].GetYaxis().SetRangeUser(ymin, ymax)
    rgraphs[0].GetXaxis().SetRangeUser(xmin, xmax)
    rgraphs[0].GetXaxis().SetLabelSize(0)
    rgraphs[0].GetYaxis().SetTitle("Ratio / Lumi")
    rgraphs[0].GetYaxis().SetTitleOffset(.75)
    rgraphs[0].GetYaxis().SetTitleSize(textsize)
    rgraphs[0].GetYaxis().SetLabelSize(textsize)
    rgraphs[0].GetYaxis().SetNdivisions(405)
    rgraphs[0].Draw("AP")

    line1 = ROOT.TLine(xmin, 1., xmax, 1)
    line1.SetLineStyle(7)
    line1.Draw("same")

    for graph in rgraphs:
        graph.Draw("P same")

    rgraphs[0].SetTitle("")
    rgraphs[0].Draw("same")


    outstring = 'ratio'
    if run_range:
        outstring += "_run{0}to{1}".format(*run_range)

    c2.SaveAs(outDir+"/"+outstring+".png")
    c2.Close()

########## Data Acquisition ##########

# --- get Z low PU
data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])

# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])


make_ratio(data_xsec, data, run_range=(299368,302029), name="2017 C/H")
make_ratio(data_xsec, data, run_range=(302030,303434), name="2017 D/H")
make_ratio(data_xsec, data, run_range=(303434,304797), name="2017 E/H")
