import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
import shutil
import pdb
import uncertainties as unc
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cms", nargs='+', help="give the CMS csv per Measurement as input", required=True)
parser.add_argument("-s", "--saveDir", default='./', type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

plotsPerFill=False
########## Data Acquisition ##########

data = pd.concat([pd.read_csv(csv, sep=',',low_memory=False) for csv in args.cms], ignore_index=True, sort=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(by=['fill','tdate_begin','tdate_end'])

# remove rows with invalid Z rates and low statistics
#data = data[np.isfinite(data['ZRate'])]
#data = data[data['z_relstat'] < 0.05]

data['tdate'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2

# data = data.loc[data['run'] > 303824]

meta = dict()

########## Plot ##########
fills = data.drop_duplicates('fill')['fill'].values

features = (
    ('ZBBeff_mc','Z barrel-barrel (corr.) efficiency', 0.8, 1.0, 'ZBBeff_stat'),
    ('ZBEeff_mc','Z barrel-endcap (corr.) efficiency', 0.8, 1.0, 'ZBEeff_stat'),
    ('ZEEeff_mc','Z endcap-endcap (corr.) efficiency', 0.8, 1.0, 'ZEEeff_stat'),
    ('ZBBeff'  ,'Z barrel-barrel efficiency', 0.8, 1.0, 'ZBBeff_stat'),
    ('ZBEeff'  ,'Z barrel-endcap efficiency', 0.8, 1.0, 'ZBEeff_stat'),
    ('ZEEeff'  ,'Z endcap-endcap efficiency', 0.8, 1.0, 'ZEEeff_stat'),
    ('HLTeffB' ,'Muon barrel HLT efficiency',0.8, 1.0, "HLTeffB_stat"),
    ('HLTeffE' ,'Muon endcap HLT efficiency',0.8, 1.0, "HLTeffE_stat"),
    ('SeleffB' ,'Muon barrel selelction efficiency',0.95, 1.0, "SeleffB_stat"),
    ('SeleffE' ,'Muon endcap selelction efficiency',0.95, 1.0, "SeleffE_stat"),
    ('TrkeffB' ,'Muon barrel track efficiency',0.975, 1.01, "TrkeffB_stat"),
    ('TrkeffE' ,'Muon endcap track efficiency',0.975, 1.01, "TrkeffE_stat"),
    ('StaeffB' ,'Muon barrel standalone efficiency',0.95, 1.01, "StaeffB_stat"),
    ('StaeffE' ,'Muon endcap standalone efficiency',0.95, 1.01, "StaeffE_stat"),
    ('zYieldBB_purity'    ,'Z barrel-barrel Purity',0.975,1.01,'zYieldBB_purity_err'),
    ('zYieldBE_purity'    ,'Z barrel-endcap Purity',0.975,1.01,'zYieldBE_purity_err'),
    ('zYieldEE_purity'    ,'Z endcap-endcap Purity',0.975,1.01,'zYieldEE_purity_err'),
    )

##### loop over Fills and produce fill specific plots
for fill in fills:
    dFill = data.loc[data['fill'] == fill]

    subDir = outDir+"/PlotsFill_"+str(fill)
    if plotsPerFill:
        if not os.path.isdir(subDir):
            os.mkdir(subDir)

    ### Efficiency ###
    for eff, name, ymin, ymax, err in features:
        # if err not in dFill.keys():
        #     dFill[err] = np.zeros(len( dFill[eff].values))

        yyUnc = dFill[eff].apply(lambda x: unc.ufloat_fromstr(x))

        yy = dFill[eff].apply(lambda x: unc.ufloat_fromstr(x).nominal_value).values
        yyErr = dFill[eff].apply(lambda x: unc.ufloat_fromstr(x).std_dev).values

        if plotsPerFill:
            graph_Zeff = ROOT.TGraphErrors(len(dFill), dFill['tdate'].values, yy, np.zeros(len(dFill['tdate'])), yyErr)
            graph_Zeff.SetName("graph_Zeff")
            graph_Zeff.SetMarkerStyle(26)
            graph_Zeff.SetMarkerColor(ROOT.kOrange+8)
            graph_Zeff.SetFillStyle(0)
            graph_Zeff.SetFillColor(0)
            graph_Zeff.SetMarkerSize(1.5)
            graph_Zeff.GetXaxis().SetTimeDisplay(1)
            graph_Zeff.SetTitle(name+", Fill "+str(fill))
            graph_Zeff.GetYaxis().SetTitle("Efficiency")
            if(eff == 'Zfpr'):
                graph_Zeff.GetYaxis().SetTitle("b/(s + b)")
            graph_Zeff.GetYaxis().SetTitleSize(0.07)
            graph_Zeff.GetYaxis().SetTitleOffset(1.1)
            graph_Zeff.GetXaxis().SetTitle("Time")
            graph_Zeff.GetXaxis().SetTitleSize(0.06)
            graph_Zeff.GetXaxis().SetTitleOffset(0.75)
            graph_Zeff.GetXaxis().SetLabelSize(0.05)
            graph_Zeff.GetYaxis().SetLabelSize(0.05)
            graph_Zeff.GetYaxis().SetRangeUser(ymin,ymax)

            c1=ROOT.TCanvas("c1","c1",1000,600)
            c1.cd(1)
            graph_Zeff.Draw("AP")
            legend=ROOT.TLegend(0.65,0.65,0.9,0.9)
            legend.AddEntry(graph_Zeff,"CMS","pe")
            #text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
            #text1.SetNDC()
            #text1.Draw()
            c1.SaveAs(subDir+"/"+str(eff)+"_"+str(fill)+".png")
            c1.Close()

        mm = np.mean(yyUnc.values)

        # sort out points with big (> 1%) statistical uncertainty
        if mm.s/mm.n > 0.01:
            continue

        if eff in meta.keys():
            meta[eff].append(mm.n)
            meta[err].append(mm.s)
        else:
            meta[eff] = [mm.n]
            meta[err] = [mm.s]

### Efficiency of all fills###

for eff, name, ymin, ymax, err in features:
    xmin = min(fills)
    xmax = max(fills)
    xmin, xmax = xmin - (xmax-xmin)*0.05, xmax + (xmax-xmin)*0.01

    c1=ROOT.TCanvas("c1","c1",800,480)
    c1.cd(1)
    c1.SetTopMargin(0.03)
    c1.SetBottomMargin(0.10)
    c1.SetRightMargin(0.01)
    c1.SetLeftMargin(0.15)

    textsize = 28./(c1.GetWh()*c1.GetAbsHNDC())

    graph_meta = ROOT.TGraphErrors(len(fills),fills.astype(float), np.array(meta[eff]), np.zeros(len(fills)), np.array(meta[err]))
    graph_meta.SetName("graph_meta")
    graph_meta.SetMarkerStyle(24)
    graph_meta.SetMarkerColor(1)#(ROOT.kOrange+8)
    graph_meta.SetFillStyle(0)
    graph_meta.SetFillColor(0)
    graph_meta.SetMarkerSize(0.5)
    graph_meta.SetTitle(name+", Fill averages ")
    graph_meta.GetYaxis().SetTitle(name)
    # if(eff == 'Zfpr'):
    #     graph_meta.GetYaxis().SetTitle("b/(s + b)")
    graph_meta.GetYaxis().SetTitleSize(textsize)
    graph_meta.GetYaxis().SetTitleOffset(1.2)
    graph_meta.GetXaxis().SetTitle("Fill")
    graph_meta.GetXaxis().SetTitleSize(textsize)
    graph_meta.GetXaxis().SetTitleOffset(0.75)
    graph_meta.GetXaxis().SetLabelSize(textsize)
    graph_meta.GetYaxis().SetLabelSize(textsize)
    graph_meta.GetYaxis().SetRangeUser(ymin,ymax)
    graph_meta.GetXaxis().SetRangeUser(xmin,xmax)

    legend=ROOT.TLegend(0.55,0.9,0.98,0.95)
    if sum(meta[err]) != 0:
        legend.AddEntry(graph_meta,"Measurement (#pm Stat.)","pe")
    else:
        legend.AddEntry(graph_meta,"Measurement","p")
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(textsize)

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(textsize)
    latex.SetTextAlign(11)


    graph_meta.Draw("AP")
    graph_meta.SetTitle("")


    latex.DrawLatex(0.26, 0.91, "Work in progress")
    latex.DrawLatex(0.18, 0.85, "138 fb^{-1} at 13 TeV")
    # latex.DrawLatex(0.18, 0.81, " ".join(name.split(" ")[:-1]))

    latex.SetTextFont(62)
    latex.DrawLatex(0.18, 0.91, 'CMS')

    legend.Draw()

    #text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text1.SetNDC()
    #text1.Draw()
    c1.SaveAs(outDir+"/allSummary_"+str(eff)+".png")
    c1.Close()
