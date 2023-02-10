###################################################################
# plot_chi2.py
# histograms of chi2 values from the fits for the Z counting measurement
###################################################################
import os,sys
import ROOT
import argparse
import pandas as pd
import numpy as np
import pdb
import matplotlib as mpl

latex = ROOT.TLatex()
latex.SetNDC()

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from python.utils import to_DateTime

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("-r","--rates", required=True, type=str, help="csv file with z rates per measurement")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
parser.add_argument("-y","--year",  default="Run-II",  type=str, help="give data taking era for labeling")
args = parser.parse_args()

year = args.year

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)


########## Data Acquisition ##########

# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(['fill','run','beginTime','endTime'])

data["beginTime"] = mpl.dates.date2num(data["beginTime"].apply(to_DateTime))
data["endTime"] = mpl.dates.date2num(data["endTime"].apply(to_DateTime))

data['time'] = (data['beginTime']+data['endTime'])/2


#sort out points with low statistics
#data = data[data['z_relstat'] < 0.05]

data = data.replace([np.inf, -np.inf], np.nan).dropna().dropna()


########## Plot ##########

### chi2 values of each category
for category in filter(lambda x: "chi2" in x,data.keys()):
    print("Now {0}".format(category))

    # only consider fits != -1 as those are dummy values
    yy = data[data[category] != -1][category]
    xx = data[data[category] != -1]['time']

    # list of fits that might have failed, based on chi2_threshold
    print("Runs that might have a failed fit:")
    failed_fits = []
    chi2_threshold = 5.0
    for x, r, m, l in data[data[category] != -1][[category,"run","measurement", "recLumi"]].values:
        if x > 5.0 or x < 0.:
            print(int(r),int(m),x, l)

    graph_chi2 = ROOT.TGraph(len(xx.values),xx.values,yy.values)
    graph_chi2.SetName("graph_chi2")
    graph_chi2.SetMarkerStyle(23)
    graph_chi2.SetMarkerColor(ROOT.kAzure-4)
    graph_chi2.SetMarkerSize(1.5)
    graph_chi2.SetTitle(category)

    graph_chi2.GetYaxis().SetTitle("#chi^{2}/ndf")
    graph_chi2.GetXaxis().SetTitle("Time")
    graph_chi2.GetXaxis().SetTimeDisplay(1)
    graph_chi2.GetXaxis().SetTimeOffset(0,"gmt")
    graph_chi2.GetXaxis().SetTitleSize(0.06)
    graph_chi2.GetYaxis().SetTitleSize(0.06)
    graph_chi2.GetXaxis().SetTitleOffset(0.72)
    graph_chi2.GetYaxis().SetTitleOffset(0.8)
    graph_chi2.GetXaxis().SetLabelSize(0.05)
    graph_chi2.GetYaxis().SetLabelSize(0.05)
    #graph_chi2.GetYaxis().SetRangeUser(-0.01,0.01)
    c3=ROOT.TCanvas("c3_"+category,"c3 "+category,1000,600)
    c3.SetGrid()

    # mean, where outlier with sigma > 1 are rejected
    avg_chi2 = np.mean(yy[abs(yy - np.mean(yy)) < np.std(yy)])

    graph_chi2.Draw("AP")

    # legend=ROOT.TLegend(0.2,0.8,0.4,0.9)
    # legend.AddEntry(graph_chi2,"Measurement","p")
    # legend.Draw("same")

    #text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text.SetNDC()
    #text.Draw()
    text2=ROOT.TLatex(0.2,0.13,"avg chi2: "+str(avg_chi2))
    text2.SetNDC()
    text2.Draw()

    c3.SaveAs(outDir+"/"+category+".png")
    c3.Close()

    if "ID" in category or "IDFail" in category:
        ndf = 120 - 5   
        #120 bins - ( 
        #   2 (signal model parameter) 
        # + 3 (background model parameters) 
        # + 1 (signal normalization) 
        # + 1 (background normalization) )
    elif "StaFail" in category:
        ndf = 120 - 7   
    elif "TrkPass" in category or  "GloPass" in category:
        ndf = 130 - 5
    elif "TrkFail" in category or  "GloFail" in category:
        ndf = 130 - 7
        #130 bins - ( 
        #   2 (signal model parameter) 
        # + 3 (background model parameters) 
        # + 1 (signal normalization) 
        # + 1 (background normalization) )
    elif "Trk" in category:
        ndf = 260 - 12
    elif "HLT2" in category or "HLT1":
        ndf = 120 - 5
    else:
        ndf = 120 - 5
        
    nBins = 50
    xmin = 0.
    xmax = 2.0*avg_chi2

    th1_chi2 = ROOT.TH1F("h1_"+category, "h1 title "+category, nBins, xmin, xmax)
    th1_chi2.GetXaxis().SetTitle("#chi^{2}/ndf")
    th1_chi2.GetYaxis().SetTitle(" ")
    th1_chi2.SetTitle(" ")
    th1_chi2.GetYaxis().SetTitleFont(42)
    th1_chi2.GetXaxis().SetTitleFont(42)
    th1_chi2.SetLineColor(1)
    th1_chi2.SetLineWidth(2)
    th1_chi2.SetFillColor(0)

    for v in data[category].values:
        th1_chi2.Fill(v)

    #add bin content for underflow/overflow bins to first/last bin
    nxbins = th1_chi2.GetNbinsX()
    th1_chi2.SetBinContent(1, th1_chi2.GetBinContent(1) + th1_chi2.GetBinContent(0))
    th1_chi2.SetBinContent(nxbins, th1_chi2.GetBinContent(nxbins+1) + th1_chi2.GetBinContent(nxbins))

    # normalize to area of one
    norm = 1./th1_chi2.Integral() / th1_chi2.GetBinWidth(1)
    th1_chi2.Scale(norm)
    # create a chi2 function to plot on top of the histogram, let ndf as free parameter
    f_chi2 = ROOT.TF1("fchi2","{0}*ROOT::Math::chisquared_pdf(x*{0},[0])".format(ndf, ndf/(xmax-xmin)), xmin, xmax)
    f_chi2.SetParameter(0, ndf)
    f_chi2.SetLineWidth(2)
    f_chi2.SetLineColor(2)
    
    # th1_chi2.Fit("fchi2")

    legend=ROOT.TLegend(0.13,0.7,0.3,0.85)
    legend.SetTextSize(0.04)
    legend.SetTextAlign(12)
    legend.SetTextFont(42)
    legend.AddEntry(f_chi2,"chi2(x,{0})".format(ndf),"l")
    legend.AddEntry(th1_chi2,"fits","f")

    canvas=ROOT.TCanvas("canvas_"+category,"canvas_"+category,1000,600)
    canvas.SetLeftMargin(0.1)
    canvas.SetRightMargin(0.03)
    canvas.SetTopMargin(0.055)
    canvas.SetBottomMargin(0.1)

    textsize = 24./(canvas.GetWh()*canvas.GetAbsHNDC())
    
    th1_chi2.SetMaximum(1.01*max(f_chi2.GetMaximum(), th1_chi2.GetMaximum()))
    
    th1_chi2.GetXaxis().SetTitleSize(textsize*1.2)
    th1_chi2.GetYaxis().SetTitleSize(textsize*1.2)
    th1_chi2.GetXaxis().SetTitleOffset(1.)
    th1_chi2.GetYaxis().SetTitleOffset(0.75)
    th1_chi2.GetXaxis().SetLabelSize(textsize*1.2)
    th1_chi2.GetYaxis().SetLabelSize(textsize*1.2)

    th1_chi2.Draw('HIST')
    f_chi2.Draw("same")
    legend.Draw("same")

    latex.SetTextFont(42)
    latex.SetTextSize(textsize*1.2)
    latex.SetTextAlign(31)

    latex.DrawLatex(0.97, 0.95, year)

    latex.SetTextAlign(11)
    latex.DrawLatex(0.5, 0.95, category.replace("_"," ").replace("chi2"," "))
    latex.DrawLatex(0.20, 0.87, "Work in progress")

    latex.SetTextFont(62)
    latex.DrawLatex(0.13, 0.87, 'CMS')

    canvas.SaveAs(outDir+"/hist_"+category+".png")
    canvas.SaveAs(outDir+"/hist_"+category+".pdf")
    canvas.Close()
