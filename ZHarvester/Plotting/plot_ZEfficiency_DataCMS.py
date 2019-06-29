import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas
import numpy as np
import shutil
import pdb

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cms", default="nothing", type=str, help="give the CMS csv as input")
parser.add_argument("-s", "--saveDir", default='./', type=str, help="give output dir")
args = parser.parse_args()

if args.cms=="nothing":
	print "please provide cms input files"
	sys.exit()

outDir = args.saveDir

########## Data Acquisition ##########

data = pandas.read_csv(str(args.cms), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(['fill','tdate_begin','tdate_end'])

# remove rows with invalid Z rates and low statistics
data = data[np.isfinite(data['ZRate'])]
data = data[data['delZCount'] > 1000.]

data['tdate'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2

meta = dict()

########## Plot ##########
fills = data.drop_duplicates('fill')['fill'].values

##### loop over Fills and produce fill specific plots
for fill in fills:
    dFill = data.loc[data['fill'] == fill]
        
    shutil.rmtree(args.saveDir+"PlotsFill_"+str(fill), ignore_errors=True)
    os.makedirs(args.saveDir+"PlotsFill_"+str(fill))
	
    ### Efficiency ###

    for eff, name, ymin, ymax in (('ZMCeff'  ,'corrected Z-Reconstruction efficitency', 0.75, 1.0),
                  ('ZMCeffBB','corrected Z-BB-Reconstruction efficitency', 0.75, 1.0),
                  ('ZMCeffBE','corrected Z-BE-Reconstruction efficitency', 0.75, 1.0),
                  ('ZMCeffEE','corrected Z-EE-Reconstruction efficitency', 0.75, 1.0),
                  ('Zeff'    ,'Z-Reconstruction efficiency', 0.75, 1.0),
                  ('ZBBeff'  ,'Z-BB-Reconstruction efficitency', 0.75, 1.0),
                  ('ZBEeff'  ,'Z-BE-Reconstruction efficitency', 0.75, 1.0),
                  ('ZEEeff'  ,'Z-EE-Reconstruction efficitency', 0.75, 1.0),
                  ('HLTeffB' ,'Muon HLT-B efficiency',0.8, 1.0 ),
                  ('HLTeffE' ,'Muon HLT-E efficiency',0.8, 1.0),
                  ('SITeffB' ,'Muon SIT-B efficiency',0.9, 1.0),
                  ('SITeffE' ,'Muon SIT-E efficiency',0.9, 1.0),
                  ('GloeffB' ,'Muon Glo-B efficiency',0.9, 1.0),
                  ('GloeffE' ,'Muon Glo-E efficiency',0.9, 1.0),
                  ('StaeffB' ,'Muon Sta-B efficiency',0.9, 1.0),
                  ('StaeffE' ,'Muon Sta-E efficiency',0.9, 1.0),
                  ('TrkeffB' ,'Muon Trk-B efficiency',0.95,1.01),
                  ('TrkeffE' ,'Muon Trk-E efficiency',0.95,1.01),
                  ('Zfpr'    ,'Z fake rate',0.0,0.1),
                 ):
        graph_Zeff = ROOT.TGraph(len(dFill),dFill['tdate'].values,dFill[eff].values )
        graph_Zeff.SetName("graph_Zeff")
        graph_Zeff.SetMarkerStyle(22)
        graph_Zeff.SetMarkerColor(ROOT.kOrange+8)
        graph_Zeff.SetFillStyle(0)
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
        c1.SaveAs(outDir+"PlotsFill_"+str(fill)+"/"+str(eff)+"_"+str(fill)+".png")
        c1.Close()
        if eff in meta.keys():
            meta[eff].append(np.mean(dFill[eff].values))
        else:
            meta[eff] = [np.mean(dFill[eff].values)]

### Efficiency of all fills###

for eff, name, ymin, ymax in (('ZMCeff'  ,'corrected Z-Reconstruction efficitency', 0.75, 1.0),
                  ('ZMCeffBB','corrected Z-BB-Reconstruction efficitency', 0.75, 1.0),
                  ('ZMCeffBE','corrected Z-BE-Reconstruction efficitency', 0.75, 1.0),
                  ('ZMCeffEE','corrected Z-EE-Reconstruction efficitency', 0.75, 1.0),
                  ('Zeff'    ,'Z-Reconstruction efficiency', 0.75, 1.0),
                  ('ZBBeff'  ,'Z-BB-Reconstruction efficitency', 0.75, 1.0),
                  ('ZBEeff'  ,'Z-BE-Reconstruction efficitency', 0.75, 1.0),
                  ('ZEEeff'  ,'Z-EE-Reconstruction efficitency', 0.75, 1.0),
                  ('HLTeffB' ,'Muon HLT-B efficiency',0.8, 1.0 ),
                  ('HLTeffE' ,'Muon HLT-E efficiency',0.8, 1.0),
                  ('SITeffB' ,'Muon SIT-B efficiency',0.9, 1.0),
                  ('SITeffE' ,'Muon SIT-E efficiency',0.9, 1.0),
                  ('GloeffB' ,'Muon Glo-B efficiency',0.9, 1.0),
                  ('GloeffE' ,'Muon Glo-E efficiency',0.9, 1.0),
                  ('StaeffB' ,'Muon Sta-B efficiency',0.9, 1.0),
                  ('StaeffE' ,'Muon Sta-E efficiency',0.9, 1.0),
                  ('TrkeffB' ,'Muon Trk-B efficiency',0.95,1.01),
                  ('TrkeffE' ,'Muon Trk-E efficiency',0.95,1.01),
                  ('Zfpr'    ,'Z fake rate',0.0,0.1),
                 ):
    graph_meta = ROOT.TGraph(len(fills),fills.astype(float),np.array(meta[eff]))
    graph_meta.SetName("graph_meta")
    graph_meta.SetMarkerStyle(22)
    graph_meta.SetMarkerColor(ROOT.kOrange+8)
    graph_meta.SetFillStyle(0)
    graph_meta.SetMarkerSize(1.5)
    graph_meta.SetTitle(name+", Fill averages ")
    graph_meta.GetYaxis().SetTitle("Efficiency")
    if(eff == 'Zfpr'):
        graph_meta.GetYaxis().SetTitle("b/(s + b)")
    graph_meta.GetYaxis().SetTitleSize(0.07)
    graph_meta.GetYaxis().SetTitleOffset(1.1)
    graph_meta.GetXaxis().SetTitle("Fill")
    graph_meta.GetXaxis().SetTitleSize(0.06)
    graph_meta.GetXaxis().SetTitleOffset(0.75)
    graph_meta.GetXaxis().SetLabelSize(0.05)
    graph_meta.GetYaxis().SetLabelSize(0.05)
    graph_meta.GetYaxis().SetRangeUser(ymin,ymax)

    c1=ROOT.TCanvas("c1","c1",1000,600)
    c1.cd(1)
    graph_meta.Draw("AP")
    legend=ROOT.TLegend(0.65,0.65,0.9,0.9)
    legend.AddEntry(graph_meta,"CMS","pe")
    #text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text1.SetNDC()
    #text1.Draw()
    c1.SaveAs(outDir+"/allSummary_"+str(eff)+".png")
    c1.Close()

