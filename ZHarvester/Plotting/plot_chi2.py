import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
import shutil

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
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()


outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)


########## Data Acquisition ##########

# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(['fill','tdate_begin','tdate_end'])

data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2

#sort out points with low statistics
#data = data[data['z_relstat'] < 0.05]

data = data.replace([np.inf, -np.inf], np.nan).dropna().dropna()


########## Plot ##########

### chi2 values of each category
for c in ('HLTeffB_chi2pass', 'HLTeffB_chi2fail', 'HLTeffE_chi2pass', 'HLTeffE_chi2fail',
          'SeleffB_chi2pass', 'SeleffB_chi2fail', 'SeleffE_chi2pass', 'SeleffE_chi2fail',
          'GloeffB_chi2pass', 'GloeffB_chi2fail', 'GloeffE_chi2pass', 'GloeffE_chi2fail',
          'zYieldBB_chi2', 'zYieldBE_chi2', 'zYieldEE_chi2'
         ):

    graph_chi2 = ROOT.TGraph(len(data),data['time'].values,data[c].values)
    graph_chi2.SetName("graph_chi2")
    graph_chi2.SetMarkerStyle(23)
    graph_chi2.SetMarkerColor(ROOT.kAzure-4)
    graph_chi2.SetMarkerSize(1.5)
    graph_chi2.SetTitle(c)

    graph_chi2.GetYaxis().SetTitle("#chi^{2}/ndf")
    graph_chi2.GetXaxis().SetTitle("Time")
    graph_chi2.GetXaxis().SetTimeDisplay(1)
    graph_chi2.GetXaxis().SetTimeOffset(0,"gmt")
    graph_chi2.GetXaxis().SetTitleSize(0.06)
    graph_chi2.GetYaxis().SetTitleSize(0.06)
    graph_chi2.GetXaxis().SetTitleOffset(0.72)
    graph_chi2.GetYaxis().SetTitleOffset(1.1)
    graph_chi2.GetXaxis().SetLabelSize(0.05)
    graph_chi2.GetYaxis().SetLabelSize(0.05)
    #graph_chi2.GetYaxis().SetRangeUser(-0.01,0.01)
    c3=ROOT.TCanvas("c3_"+c,"c3 "+c,1000,600)
    c3.SetGrid()

    # mean, where outlier with sigma > 1 are rejected
    avg_chi2 = np.mean(data[c][abs(data[c] - np.mean(data[c])) < np.std(data[c])])

    graph_chi2.Draw("AP")

    legend=ROOT.TLegend(0.2,0.8,0.4,0.9)
    legend.AddEntry(graph_chi2,"Measurement","p")
    legend.Draw("same")

    #text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text.SetNDC()
    #text.Draw()
    text2=ROOT.TLatex(0.2,0.23,"avg chi2: "+str(avg_chi2))
    text2.SetNDC()
    text2.Draw()

    c3.SaveAs(outDir+"/"+c+".png")
    c3.Close()

    th1_chi2 = ROOT.TH1F("h1_"+c, "h1 title "+c, 20, 0., 2.0*avg_chi2)
    th1_chi2.GetXaxis().SetTitle("#chi^{2}/ndf")
    th1_chi2.GetYaxis().SetTitle("fits")
    th1_chi2.GetXaxis().SetTitleSize(0.06)
    th1_chi2.GetYaxis().SetTitleSize(0.06)
    th1_chi2.GetXaxis().SetTitleOffset(0.72)
    th1_chi2.GetYaxis().SetTitleOffset(1.1)
    th1_chi2.GetXaxis().SetLabelSize(0.05)
    th1_chi2.GetYaxis().SetLabelSize(0.05)
    th1_chi2.SetTitle(c)
    th1_chi2.SetLineColor(1)
    th1_chi2.SetLineWidth(2)
    th1_chi2.SetFillColor(0)

    for v in data[c].values:
        th1_chi2.Fill(v)

    #f_chi2 = ROOT.TF1("fchi2","{0}*TMath::Prob(x,1)".format(th1_chi2.Integral()*(2.0*avg_chi2)/20.),0.0,2.0*avg_chi2)
    #f_chi2.SetLineWidth(2)
    #f_chi2.SetLineColor(2)

    #legend=ROOT.TLegend(0.8,0.8,0.9,0.9)
    #legend.AddEntry(f_chi2,"chi2(x,1)","l")
    #legend.AddEntry(th1_chi2,"fits","f")

    c4=ROOT.TCanvas("c4_"+c,"c4_"+c,1000,600)

    th1_chi2.Draw('HIST')
    #f_chi2.Draw("same")
    #legend.Draw("same")


    c4.SaveAs(outDir+"/hist_"+c+".png")
    c4.Close()
