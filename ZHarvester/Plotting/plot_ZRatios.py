import os,sys
import ROOT
import argparse
import pandas as pd
import numpy as np

import pdb

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("--rates1", required=True, type=str, help="Nominator csv file with z rates perLS")
parser.add_argument("--rates2",  required=True, type=str, help="Denominator second csv file with z rates perLS")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)


def make_ratio(dataNom, dataDenom, run_range_Nom=None, run_range_Denom=None,
    lumiUnc=0.013, name="", lumi_name='lumiRec'):

    if run_range_Nom:
        dataNom = dataNom.query("run >= {0} & run <= {1}".format(*run_range_Nom))

    if run_range_Denom:
        dataDenom = dataDenom.query("run >= {0} & run <= {1}".format(*run_range_Denom))

    rLumi = dataNom[lumi_name].sum() / dataDenom[lumi_name].sum()

    rZ_BB = dataNom['zDelBB_mc'].sum() / dataDenom['zDelBB_mc'].sum()
    rZ_BE = dataNom['zDelBE_mc'].sum() / dataDenom['zDelBE_mc'].sum()
    rZ_EE = dataNom['zDelEE_mc'].sum() / dataDenom['zDelEE_mc'].sum()
    rZ_tot = dataNom['zDel_mc'].sum() / dataDenom['zDel_mc'].sum()

    # uncertainty on pileup correction: assume 100% uncertainty
    rZ_BB_err_PU = abs(rZ_BB - dataNom['zDelBB'].sum() / dataDenom['zDelBB'].sum())
    rZ_BE_err_PU = abs(rZ_BE - dataNom['zDelBE'].sum() / dataDenom['zDelBE'].sum())
    rZ_EE_err_PU = abs(rZ_EE - dataNom['zDelEE'].sum() / dataDenom['zDelEE'].sum())
    rZ_tot_err_PU = abs(rZ_tot - dataNom['zDel'].sum() / dataDenom['zDel'].sum())

    rZ_BB_err = rZ_BB * 1. / np.sqrt(dataDenom['zYieldBB'].sum())
    rZ_BE_err = rZ_BE * 1. / np.sqrt(dataDenom['zYieldBE'].sum())
    rZ_EE_err = rZ_EE * 1. / np.sqrt(dataDenom['zYieldEE'].sum())
    rZ_tot_err = rZ_tot * 1. / (dataDenom['zDelBB_mc'].sum() + dataDenom['zDelBE_mc'].sum() + dataDenom['zDelEE_mc'].sum()) * np.sqrt(dataDenom['zYieldBB'].sum() + dataDenom['zYieldBE'].sum() + dataDenom['zYieldEE'].sum())

    points = np.array([rLumi, rZ_BB, rZ_BE, rZ_EE, rZ_tot])
    points_err = np.array([rLumi*lumiUnc, rZ_BB_err, rZ_BE_err, rZ_EE_err, rZ_tot_err])
    points_err2 = np.array([0., rZ_BB_err_PU, rZ_BE_err_PU, rZ_EE_err_PU, rZ_tot_err_PU])

    ########## Plot ##########
    xmin = 0.9
    xmax = 1.5

    graphs = []
    graphs2 = []
    for i, (ipoint, ierr, ierr2, nam, ptr) in enumerate(
        zip(points, points_err, points_err2,
            ['Lumi','Z BB','Z BE','Z EE','Z total'],
            [20, 21, 22, 23, 34],
        )):
        graph=ROOT.TGraphErrors(1, np.array([0.1*i+1]), np.array([ipoint]),  np.array([0.]), np.array(np.sqrt(np.array(ierr)**2 + np.array(ierr2)**2)))
        graph2=ROOT.TGraphErrors(1, np.array([0.1*i+1]), np.array([ipoint]),  np.array([0.]), np.array(ierr2))
        graph.SetName(nam)
        graph.SetTitle(nam)
        graph.SetMarkerStyle(ptr)
        graph2.SetMarkerStyle(ptr)
        graph.SetMarkerColor(i+1)
        graph2.SetMarkerColor(i+1)
        graph.SetFillStyle(1001)
        graph2.SetFillStyle(1001)
        graph.SetMarkerSize(1.5)
        graph2.SetMarkerSize(1.5)
        graph.SetLineColor(1)
        graph2.SetLineColor(2)
        graphs.append(graph)
        graphs2.append(graph2)

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
    ymax = ymin + 1.25 * (max(points) - ymin)

    graphs[0].GetYaxis().SetRangeUser(ymin, ymax)
    graphs[0].GetXaxis().SetRangeUser(xmin, xmax)
    graphs[0].GetXaxis().SetLabelSize(0)
    graphs[0].GetYaxis().SetTitle("Ratio")
    graphs[0].GetYaxis().SetTitleOffset(1.4)
    graphs[0].GetYaxis().SetTitleSize(textsize)
    graphs[0].GetYaxis().SetLabelSize(textsize)
    graphs[0].Draw("AP")

    legend=ROOT.TLegend(0.75,0.55,0.98,0.85)

    latex.SetTextSize(textsize)
    latex.SetTextFont(42)

    latex.DrawLatex(0.2, 0.91, name)

    latex.SetTextAlign(11)
    latex.DrawLatex(0.35, 0.81, "Preliminary")

    latex.SetTextAlign(11)
    latex.SetTextFont(62)
    latex.DrawLatex(0.23, 0.81, 'CMS')

    for graph in graphs:
        graph.Draw("P same")
        legend.AddEntry(graph,"","pe")

    for graph in graphs2:
        graph.Draw("E same")

    legend.SetTextFont(42)
    legend.SetTextSize(textsize)

    legend.Draw("same")

    graphs[0].SetTitle("")
    graphs[0].Draw("same")

    ### ratio ###

    points_err /= points[0]
    points_err2 /= points[0]
    points /= points[0]

    rgraphs = []
    rgraphs2 = []
    for i, (ipoint, ierr, ierr2, nam, ptr) in enumerate(
        zip(points, points_err, points_err2,
            ['Lumi','Z BB','Z BE','Z EE','Z total'],
            [20, 21, 22, 23, 34],
        )):
        graph=ROOT.TGraphErrors(1, np.array([0.1*i+1]), np.array([ipoint]),  np.array([0.]), np.array(np.sqrt(np.array(ierr)**2 + np.array(ierr2)**2)))
        graph2=ROOT.TGraphErrors(1, np.array([0.1*i+1]), np.array([ipoint]),  np.array([0.]), np.array(ierr2))
        graph.SetName(nam)
        graph.SetTitle(nam)
        graph.SetMarkerStyle(ptr)
        graph2.SetMarkerStyle(ptr)
        graph.SetMarkerColor(i+1)
        graph2.SetMarkerColor(i+1)
        graph.SetFillStyle(1001)
        graph2.SetFillStyle(1001)
        graph.SetMarkerSize(1.5)
        graph2.SetMarkerSize(1.5)
        graph.SetLineColor(1)
        graph2.SetLineColor(2)
        rgraphs.append(graph)
        rgraphs2.append(graph2)

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

    for graph in rgraphs2:
        graph.Draw("E same")

    rgraphs[0].SetTitle("")
    rgraphs[0].Draw("same")

    outstring = 'ratio'
    if run_range_Nom:
        outstring += "_run{0}to{1}".format(*run_range_Nom)
    if run_range_Denom:
        outstring += "_run{0}to{1}".format(*run_range_Denom)

    c2.SaveAs(outDir+"/"+outstring+".png")
    c2.Close()

########## Data Acquisition ##########

# --- z luminosity
dataNom = pd.read_csv(str(args.rates1), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])

# --- get Z low PU
dataDenom = pd.read_csv(str(args.rates2), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])

dataNom['zDel_mc'] = dataNom['zDelBB_mc'] + dataNom['zDelBE_mc'] + dataNom['zDelEE_mc']
dataNom['zDel'] = dataNom['zDelBB'] + dataNom['zDelBE'] + dataNom['zDelEE']
dataDenom['zDel_mc'] = dataDenom['zDelBB_mc'] + dataDenom['zDelBE_mc'] + dataDenom['zDelEE_mc']
dataDenom['zDel'] = dataDenom['zDelBB'] + dataDenom['zDelBE'] + dataDenom['zDelEE']

# pdb.set_trace()
# sort out lumi section withou any counts
dataNom = dataNom.query('zDel_mc != 0')
dataDenom = dataDenom.query('zDel_mc != 0')


# make_ratio(dataDenom, dataNom, run_range=(297046,299329), name="2017 B/H")
# make_ratio(dataDenom, dataNom, run_range=(299368,302029), name="2017 C/H")
# make_ratio(dataDenom, dataNom, run_range=(302030,303434), name="2017 D/H")
# make_ratio(dataDenom, dataNom, run_range=(303434,304797), name="2017 E/H")
# make_ratio(dataDenom, dataNom, run_range=(305040,306462), name="2017 F/H")

#make_ratio(dataDenom, dataDenom, run_range_Denom=(317080,319310), run_range_Nom=(315252,316995), name="2018 A / 2018 B", lumiUnc=0.)

make_ratio(dataNom, dataDenom,
    run_range_Nom=(315252,320065),
    run_range_Denom=(297046,306462),
    name="2018 ABC / 2017 B-F",
    lumiUnc=np.sqrt(0.022**2 + 0.015**2),
    lumi_name='recorded(/pb)')

# make_ratio(dataNom, dataDenom,
#     name="2018 ABC / 2017 H",
#     run_range_Nom=(315252,320065),
#     lumiUnc=np.sqrt(0.015**2 + 0.015**2),
#     lumi_name='recorded(/pb)')

# make_ratio(dataNom, dataDenom,
#     run_range_Nom=(297046,306462),
#     name="2017 B-F / 2017 H",
#     lumiUnc=0.013,
#     lumi_name='recorded(/pb)')

# make_ratio(dataNom, dataDenom,
#     run_range_Nom=(303434,306462),
#     run_range_Denom=(303434,306462),
#     name="2017 EF(0.1) / 2017 EF(0.2)",
#     lumiUnc=0.013,
#     lumi_name='recorded(/pb)')
