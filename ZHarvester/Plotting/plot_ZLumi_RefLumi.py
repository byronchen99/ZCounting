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
parser.add_argument("--refLumi",  required=True, type=str, help="give a ByLs.csv as input for reference Luminosity")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Constants ##########

secPerLS=float(23.3)


### For plot labeling ###

currentYear = 2017

refLumiSource='PHYSICS'

########## Data Acquisition ##########

# --- reference luminosity
# convert reference data from ByLs.csv file
lumiFile=open(str(args.refLumi))
lumiLines=lumiFile.readlines()
data_ref = pd.read_csv(str(args.refLumi), sep=',',low_memory=False, skiprows=lambda x: lumiLines[x].startswith('#') and not lumiLines[x].startswith('#run'))
if 'recorded(/ub)' in data_ref.columns.tolist():      #convert to /pb
    data_ref['recorded(/ub)'] = data_ref['recorded(/ub)'].apply(lambda x:x / 1000000.)
    data_ref = data_ref.rename(index=str, columns={'recorded(/ub)':'recorded(/pb)' })
elif 'recorded(/fb)' in data_ref.columns.tolist():      #convert to /pb
    data_ref['recorded(/fb)'] = data_ref['recorded(/fb)'].apply(lambda x:x * 1000.)
    data_ref = data_ref.rename(index=str, columns={'recorded(/fb)':'recorded(/pb)' })

data_ref['fill'] = pd.to_numeric(data_ref['#run:fill'].str.split(':',expand=True)[1])
data_ref['run'] = pd.to_numeric(data_ref['#run:fill'].str.split(':',expand=True)[0])
data_ref['ls'] = pd.to_numeric(data_ref['ls'].str.split(':',expand=True)[0])
data_ref = data_ref.drop(['#run:fill','hltpath','source'],axis=1)
data_ref = data_ref.sort_values(['fill','run','ls','recorded(/pb)'])
data_ref = data_ref.drop_duplicates(['fill','run','ls'])

data_ref['time'] = data_ref['time'].apply(lambda x: to_RootTime(x,currentYear)).astype(float)

data_ref['timeE'] = data_ref['time'] - data_ref['time']
data_ref['dLRec(/pb)'] = data_ref['recorded(/pb)']/secPerLS
data_ref = data_ref[['fill','dLRec(/pb)','time']]		#Keep only what you need

# --- get Z xsec
data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
xsecBB = sum(data_xsec['zDelBB_mc'])/sum(data_xsec['lumiRec'])
xsecBE = sum(data_xsec['zDelBE_mc'])/sum(data_xsec['lumiRec'])
xsecEE = sum(data_xsec['zDelEE_mc'])/sum(data_xsec['lumiRec'])

# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(['fill','tdate_begin','tdate_end'])

data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2
data['timeE'] = (data['tdate_end'] - data['tdate_begin'])/2

data['zLumiInstBB_mc'] = data['zDelBB_mc'] / (data['timewindow'] * xsecBB)
data['zLumiInstBE_mc'] = data['zDelBE_mc'] / (data['timewindow'] * xsecBE)
data['zLumiInstEE_mc'] = data['zDelEE_mc'] / (data['timewindow'] * xsecEE)
data['zLumiInst_mc'] = (data['zDelBB_mc'] + data['zDelBE_mc'] + data['zDelEE_mc']) / (data['timewindow'] * (xsecBB + xsecBE + xsecEE))

data['zLumiInstBB_mc_stat'] = data['zLumiInstBB_mc'] * data['zBB_relstat']
data['zLumiInstBE_mc_stat'] = data['zLumiInstBE_mc'] * data['zBE_relstat']
data['zLumiInstEE_mc_stat'] = data['zLumiInstEE_mc'] * data['zEE_relstat']
# BB and EE uncorrelated, but partially correlated with BE
data['zLumiInst_mc_stat'] = data['zLumiInst_mc'] * (np.sqrt(data['zEE_relstat']**2 + data['zBB_relstat']**2) + data['zBE_relstat'])


data = data.replace([np.inf, -np.inf], np.nan).dropna().dropna()


########## Plot ##########

for suffix, zLumi, zLumiErr in (
    ("inclusive", 'zLumiInst_mc',   'zLumiInst_mc_stat'),
    ("BB",        'zLumiInstBB_mc', 'zLumiInstBB_mc_stat'),
    ("BE",        'zLumiInstBE_mc', 'zLumiInstBE_mc_stat'),
    ("EE",        'zLumiInstEE_mc', 'zLumiInstEE_mc_stat'),
):
    ##### loop over Fills and produce fill specific plots
    for fill in data.drop_duplicates('fill')['fill'].values:
        dFill = data.loc[data['fill'] == fill]

        # startTime=dFill.iloc[0]['time']
        # endTime = dFill.iloc[-1]['time']

        if not os.path.isdir(outDir):
            os.mkdir(outDir+"/PlotsFill_"+str(fill))

        ### Z Luminosity vs reference Luminosity ###

        dFill_ref = data_ref.loc[data_ref['fill'] == fill]
        graph_Lumi=ROOT.TGraph(len(dFill_ref),dFill_ref['time'].values,  dFill_ref['dLRec(/pb)'].values)
        graph_Lumi.SetName("graph_Lumi")
        graph_Lumi.SetMarkerStyle(23)
        graph_Lumi.SetMarkerColor(ROOT.kBlue+8)
        graph_Lumi.SetFillStyle(0)
        graph_Lumi.SetMarkerSize(1.5)

        graph_ZLumi=ROOT.TGraphErrors(len(dFill),dFill['time'].values,dFill[zLumi].values, dFill['timeE'].values, dFill[zLumiErr].values )
        graph_ZLumi.SetName("graph_ZLumi")
        graph_ZLumi.SetMarkerStyle(22)
        graph_ZLumi.SetMarkerColor(ROOT.kOrange+8)
        graph_ZLumi.SetFillStyle(0)
        graph_ZLumi.SetMarkerSize(1.5)

        graph_ZLumi.GetXaxis().SetTimeDisplay(1)
        graph_ZLumi.GetYaxis().SetTitle("inst. Luminosity [1/pb]")
        graph_ZLumi.GetYaxis().SetTitleSize(0.07)
        graph_ZLumi.GetYaxis().SetTitleOffset(1.1)
        graph_ZLumi.GetXaxis().SetTitle("Time")
        graph_ZLumi.GetXaxis().SetTitleSize(0.06)
        graph_ZLumi.GetXaxis().SetTitleOffset(0.75)
        graph_ZLumi.GetXaxis().SetLabelSize(0.05)
        graph_ZLumi.GetYaxis().SetLabelSize(0.05)
        #graph_ZLumi.GetXaxis().SetRangeUser(startTime,endTime)

        graph_ZLumi.SetTitle(suffix+" inst. Luminosity, Fill "+str(fill))

        c2=ROOT.TCanvas("c2","c2",1000,600)
        c2.cd(1)
        graph_ZLumi.Draw("AP")
        graph_Lumi.Draw("l same")

        legend=ROOT.TLegend(0.65,0.65,0.9,0.8)
        legend.AddEntry(graph_ZLumi,"Z Lumi","pe")
        legend.AddEntry(graph_Lumi,refLumiSource+" Lumi","l")
        legend.Draw("same")

        graph_ZLumi.Draw("P same")

        #text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        #text1.SetNDC()
        #text1.Draw()
        c2.SaveAs(outDir+"/PlotsFill_"+str(fill)+"/ZLumi"+str(fill)+"_"+suffix+".png")
        c2.Close()
