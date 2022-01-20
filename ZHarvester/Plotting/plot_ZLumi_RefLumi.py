import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
import uncertainties as unc
import pdb

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime, getMCCorrection, cms, preliminary, text, workinprogress
from python.corrections import apply_muon_prefire, apply_ECAL_prefire

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

ROOT.gStyle.SetErrorX(0.)
ROOT.gStyle.SetTextFont(42)
ROOT.gStyle.SetTitleFont(42)
ROOT.gStyle.SetLabelFont(42)
ROOT.gStyle.SetLegendFont(42)
# ROOT.gStyle.SetLegendTextSize(0.6)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetEndErrorSize(2)

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--rates", required=True, type=str, help="csv file with z rates per measurement")
parser.add_argument("-x", "--xsec",  default=None, type=str,
    help="csv file with z rates per measurement where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("-l", "--refLumi",  required=True, type=str, help="give a ByLs.csv as input for reference Luminosity")
parser.add_argument("-s", "--saveDir",  default='./',  type=str, help="give output dir")
parser.add_argument("-y", "--year",  default=2017, type=int, help="give a year for calculation of time")
parser.add_argument("-f", "--fill", nargs="*",  type=int,  default=[], help="specify a single fill to plot")
parser.add_argument("-m", "--mcCorrections",  default=None, type=str, help="specify MC correction file")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Constants ##########

secPerLS=float(23.3)

# normalize Z luminosity to PHYSICS luminosity for each fill, if all regions are shown, normalized against inclusive
normalized = True
margin_left   = 0.125
margin_right  = 0.02
margin_top    = 0.02
margin_bottom = 0.35
### For plot labeling ###

currentYear = args.year

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

if args.fill != []:
    data_ref = data_ref.loc[data_ref['fill'].isin(data_ref.fill)]

data_ref['time'] = data_ref['time'].apply(lambda x: to_RootTime(x,currentYear)).astype(float)

data_ref['timeE'] = data_ref['time'] - data_ref['time']
data_ref['dLRec(/pb)'] = data_ref['recorded(/pb)']/secPerLS
data_ref = data_ref[['fill','dLRec(/pb)','time', 'recorded(/pb)', 'avgpu']]		#Keep only what you need

# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False) #, skiprows=[1,2,3,4,5])
if args.fill != []:
    data = data.loc[data['fill'].isin(args.fill)]

data = data.sort_values(['fill','tdate_begin','tdate_end'])

data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2
data['timeE'] = (data['tdate_end'] - data['tdate_begin'])/2

def unorm(x):
    # for counting experiments: define ufloat with poisson uncertainty
    return unc.ufloat(x, np.sqrt(abs(x)))

# data['zDelBB_mc'] = data['zDelBB'].apply(lambda x: unc.ufloat_fromstr(x))
# data['zDelBE_mc'] = data['zDelBE'].apply(lambda x: unc.ufloat_fromstr(x))
# data['zDelEE_mc'] = data['zDelEE'].apply(lambda x: unc.ufloat_fromstr(x))

# take uncertainties from uncorrected zDel
data['zDelBB_mc'] = data['zDelBB'].apply(lambda x: unc.ufloat_fromstr(x)/unc.ufloat_fromstr(x).n) * data['zDelBB_mc']
data['zDelBE_mc'] = data['zDelBE'].apply(lambda x: unc.ufloat_fromstr(x)/unc.ufloat_fromstr(x).n) * data['zDelBE_mc']
data['zDelEE_mc'] = data['zDelEE'].apply(lambda x: unc.ufloat_fromstr(x)/unc.ufloat_fromstr(x).n) * data['zDelEE_mc']


# --- get Z xsec
if args.xsec:
    data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])

    data_xsec['zDelBB_mc'] = data_xsec['zDelBB_mc'].apply(lambda x: unc.ufloat_fromstr(x))
    data_xsec['zDelBE_mc'] = data_xsec['zDelBE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
    data_xsec['zDelEE_mc'] = data_xsec['zDelEE_mc'].apply(lambda x: unc.ufloat_fromstr(x))

    xsecBB = sum(data_xsec['zDelBB_mc']) / sum(data_xsec['lumiRec'])
    xsecBE = sum(data_xsec['zDelBE_mc']) / sum(data_xsec['lumiRec'])
    xsecEE = sum(data_xsec['zDelEE_mc']) / sum(data_xsec['lumiRec'])

else:
    xsecBB = 1.
    xsecBE = 1.
    xsecEE = 1.

# --->>> prefire corrections
apply_muon_prefire(data)
apply_ECAL_prefire(data)

data['zLumi_mc'] = (data['zDelBB_mc'] + data['zDelBE_mc'] + data['zDelEE_mc']) / (xsecBB+xsecBE+xsecEE)
data['zLumiBB_mc'] = data['zDelBB_mc'] / xsecBB
data['zLumiBE_mc'] = data['zDelBE_mc'] / xsecBE
data['zLumiEE_mc'] = data['zDelEE_mc'] / xsecEE

data['zLumiInst_mc'] = data['zLumi_mc'] / data['timewindow']
data['zLumiInstBB_mc'] = data['zLumiBB_mc'] / data['timewindow']
data['zLumiInstBE_mc'] = data['zLumiBE_mc'] / data['timewindow']
data['zLumiInstEE_mc'] = data['zLumiEE_mc'] / data['timewindow']


# seconds to hours
data['time'] = data['time'] / 3600.
data['timeE'] = data['timeE'] / 3600.
data['tdate_begin'] = data['tdate_begin'] / 3600.
data['tdate_end'] = data['tdate_end'] / 3600.

data_ref['time'] = data_ref['time'] / 3600.

data = data.replace([np.inf, -np.inf], np.nan).dropna().dropna()

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'same') / w

########## Plot ##########

for suffix, zLumi, others in (
    ("inclusive", 'zLumiInst_mc',   [], ), #[("zLumiInst",   24, 2, "w/o MC"),] ),
    # ("BB",        'zLumiInstBB_mc', [("zLumiInstBB", 24, 2, "w/o MC"),] ),
    # ("BE",        'zLumiInstBE_mc', [("zLumiInstBE", 24, 2, "w/o MC"),] ),
    # ("EE",        'zLumiInstEE_mc', [("zLumiInstEE", 24, 2, "w/o MC"),] ),
    # ("All",       'All', [ ("zLumiInstBB_mc", 24, 2, "BB"), ("zLumiInstBE_mc", 25, 6, "BE"), ("zLumiInstEE_mc", 32, 4, "EE")]),
):
    ##### loop over Fills and produce fill specific plots
    for fill in data.drop_duplicates('fill')['fill'].values:
        dFill = data.loc[data['fill'] == fill]

        ### Z Luminosity vs reference Luminosity ###
        dFill_ref = data_ref.loc[data_ref['fill'] == fill]

        subdir = outDir+"/PlotsFill_"+str(fill)
        if not os.path.isdir(subdir):
            os.mkdir(subdir)

        # sum up reference lumi for each measurement (needed for ratio)
        dFill_refs = []
        for tBegin, tEnd, time, pileup in dFill[['tdate_begin','tdate_end', 'time', 'pileUp']].values:
            dFill_iRef = dFill_ref.loc[(dFill_ref['time'] >= tBegin) & (dFill_ref['time'] < tEnd)].mean()
            dFill_iRef['time'] = time
            dFill_iRef['pileUp'] = pileup
            dFill_refs.append(dFill_iRef)

        dFill_cRef = pd.concat(dFill_refs, axis=1).T

        dFill = pd.merge(dFill, dFill_cRef, "outer")

        starttime = dFill['time'].values[0]
        dFill_time = dFill['time'].values - starttime
        dFill_timeE = dFill['timeE'].values

        dFill_PU = dFill['pileUp'].values
        dFill_PUE = np.zeros(len(dFill_PU))

        yyDelL = dFill_ref['dLRec(/pb)'].values


        xxRefTime = dFill_ref['time'].values - starttime
        xxRefPileup = dFill_ref['avgpu'].values

        for xxRef, yy, xxZ, xxZErr, xLabel, xName, legend_position in (
            (xxRefTime,   yyDelL, dFill_time, dFill_timeE, "LHC Runtime [h]", "time", 1),
            (xxRefPileup, yyDelL, dFill_PU,   dFill_PUE,   "<Number of pileup interactions>", "pileup", 2)
            ):

            # sort Reference lumi
            zipped = zip(xxRef, yy)
            zipped.sort()
            xxRef, yy = zip(*zipped)

            # smoothening
            yy = moving_average(yy, 40)


            graph_Lumi=ROOT.TGraph(len(xxRef), np.array(xxRef), np.array(yy) * 1000)
            graph_Lumi.SetName("graph_Lumi")
            # graph_Lumi.SetMarkerStyle(23)
            # graph_Lumi.SetMarkerColor(ROOT.kBlue+8)
            graph_Lumi.SetFillStyle(0)
            graph_Lumi.SetLineWidth(2)
            graph_Lumi.SetLineColor(1)
            graph_Lumi.SetMarkerSize(1.5)
            graph_Lumi.SetTitle("")

            if normalized:
                norm = sum(dFill_ref['recorded(/pb)'])/sum(dFill[zLumi.replace("Inst","")])
            else:
                norm = 1.0

            graphs = []
            if others:

                for zL, ms, mc, label in others:

                    if normalized:
                        dFill_zL = dFill[zL] * sum(dFill_ref['recorded(/pb)'])/sum(dFill[zL.replace("Inst","")])
                    else:
                        dFill_zL = dFill[zL]

                    graph = ROOT.TGraph(len(dFill),
                        xxZ, dFill_zL.apply(lambda x: x.n).values * 1000)
                        # , dFill_timeE, dFill_zL.apply(lambda x: x.s).values * 1000 )

                    graph.SetMarkerStyle(ms)
                    graph.SetMarkerColor(mc)
                    graph.SetFillStyle(0)
                    graph.SetFillColor(0)
                    graph.SetMarkerSize(1.5)
                    # graph.GetXaxis().SetTimeDisplay(1)
                    graph.SetTitle("Z Lumi ({0})".format(label))

                    graphs.append(graph)

            dFill_zLumi = dFill[zLumi] * norm

            graph_ZLumi=ROOT.TGraphErrors(len(dFill),
                xxZ, dFill_zLumi.apply(lambda x: x.n).values * 1000,
                xxZErr, dFill_zLumi.apply(lambda x: x.s).values * 1000 )
            graph_ZLumi.SetName("graph_ZLumi")
            graph_ZLumi.SetTitle("")
            graph_ZLumi.SetMarkerStyle(20)
            graph_ZLumi.SetMarkerColor(2)
            graph_ZLumi.SetFillStyle(0)
            graph_ZLumi.SetFillColor(0)
            graph_ZLumi.SetMarkerSize(1)


            xmin = min((xxZ - xxZErr))
            xmax = max((xxZ + xxZErr))
            xWidth = xmax - xmin
            xmin = xmin - (xWidth * 0.05)
            xmax = xmax + (xWidth * 0.05)

            canvas=ROOT.TCanvas("canvas","canvas",800,640)
            pad1 = ROOT.TPad("pad1", "pad1", 0., 0.35, 1, 1.0)
            pad1.SetBottomMargin(0.01)
            canvas.SetTicks()
            pad1.SetLeftMargin(margin_left)
            pad1.SetRightMargin(margin_right)
            pad1.SetTopMargin(margin_top)
            pad1.SetTickx()
            pad1.SetTicky()
            pad1.Draw()
            pad1.cd()

            textsize1 = 28./(pad1.GetWh()*pad1.GetAbsHNDC())

            graph_ZLumi.GetYaxis().SetTitle("inst. Luminosity [1/nb]")
            graph_ZLumi.GetYaxis().SetTitleSize(textsize1)
            graph_ZLumi.GetYaxis().SetTitleOffset(0.9)
            graph_ZLumi.GetYaxis().SetLabelSize(textsize1)

            # graph_ZLumi.GetXaxis().SetTimeDisplay(1)
            graph_ZLumi.GetXaxis().SetTitleSize(0)
            graph_ZLumi.GetXaxis().SetLabelSize(0)
            graph_ZLumi.GetXaxis().SetRangeUser(xmin,xmax)

            graph_ZLumi.Draw("AP")

            graph_Lumi.Draw("l same")

            if legend_position == 1:
                legend=ROOT.TLegend(0.62,0.72,0.92,0.92)
            elif legend_position == 2:
                legend=ROOT.TLegend(0.62,0.08,0.92,0.32)

            legend.SetTextSize(textsize1)
            legend.AddEntry(graph_ZLumi,"Z Lumi","pe")

            for graph in graphs:
                legend.AddEntry(graph,graph.GetTitle(),"p")
                graph.SetTitle("")
                graph.Draw("P same")

            legend.AddEntry(graph_Lumi,refLumiSource+" Lumi","l")
            legend.Draw("same")

            # cms(x=0.5, y=0.88, textsize=textsize1)
            # preliminary(x=0.35, y=0.88, textsize=textsize1)
            if legend_position == 1:
                cms(x=0.28, y=0.16, textsize=textsize1)
                workinprogress(x=0.38, y=0.16, textsize=textsize1)
                text("Fill "+str(fill), x=0.28, y=0.08, textsize=textsize1)
            elif legend_position == 2:
                cms(x=0.18, y=0.88, textsize=textsize1)
                workinprogress(x=0.28, y=0.88, textsize=textsize1)
                text("Fill "+str(fill), x=0.18, y=0.80, textsize=textsize1)

            graph_ZLumi.Draw("P same")


            # --- Ratio plot
            ymin = 0.901
            ymax = 1.099

            canvas.cd()
            pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.35)
            pad2.SetLeftMargin(margin_left)
            pad2.SetRightMargin(margin_right)
            pad2.SetTopMargin(0.01)
            pad2.SetBottomMargin(margin_bottom)
            pad2.SetTickx()
            pad2.SetTicky()
            pad2.Draw("ALPF")
            pad2.cd()

            textsize2 = 28./(pad2.GetWh()*pad2.GetAbsHNDC())

            graphs_r = []
            if others:

                for zL, ms, mc, label in others:
                    if normalized:
                        dFill_zL = dFill[zL] * sum(dFill_ref['recorded(/pb)'])/sum(dFill[zL.replace("Inst","")]) / dFill['dLRec(/pb)']
                    else:
                        dFill_zL = dFill[zL] / dFill['dLRec(/pb)']

                    graph = ROOT.TGraph(len(dFill), xxZ, dFill_zL.apply(lambda x: x.n).values)
                        # , dFill_timeE, dFill_zL.apply(lambda x: x.s).values )

                    graph.SetMarkerStyle(ms)
                    graph.SetMarkerColor(1)
                    graph.SetFillStyle(0)
                    graph.SetFillColor(0)
                    graph.SetMarkerSize(1)
                    # graph.GetXaxis().SetTimeDisplay(1)
                    graph.SetTitle("Z Lumi ({0})".format(label))

                    graphs_r.append(graph)

            dFill_zLumi = dFill_zLumi / dFill['dLRec(/pb)']

            graph_rZLumi=ROOT.TGraphErrors(len(dFill),
                xxZ, dFill_zLumi.apply(lambda x: x.n).values,
                xxZErr, dFill_zLumi.apply(lambda x: x.s).values )
            graph_rZLumi.SetName("graph_ZLumi")
            graph_rZLumi.SetTitle("")
            graph_rZLumi.SetMarkerStyle(20)
            graph_rZLumi.SetMarkerColor(2)
            graph_rZLumi.SetFillStyle(0)
            graph_rZLumi.SetMarkerSize(1)

            graph_rZLumi.GetYaxis().SetNdivisions(405)
            graph_rZLumi.GetYaxis().SetTitle("Ratio")
            graph_rZLumi.GetYaxis().SetTitleOffset(0.5)
            graph_rZLumi.GetYaxis().SetTitleSize(textsize2)
            graph_rZLumi.GetYaxis().SetLabelSize(textsize2)
            graph_rZLumi.GetYaxis().CenterTitle(True)

            # graph_rZLumi.GetXaxis().SetTimeDisplay(1)
            graph_rZLumi.GetXaxis().SetTitle(xLabel)
            graph_rZLumi.GetXaxis().SetTitleSize(textsize2)
            graph_rZLumi.GetXaxis().SetTitleOffset(1.0)
            graph_rZLumi.GetXaxis().SetLabelSize(textsize2)
            graph_rZLumi.GetXaxis().SetRangeUser(xmin,xmax)

            graph_rZLumi.SetMinimum(ymin)
            graph_rZLumi.SetMaximum(ymax)
            graph_rZLumi.GetXaxis().SetRangeUser(xmin,xmax)

            graph_rZLumi.Draw("AP")


            line0 = ROOT.TLine(xmin, 1, xmax, 1)
            line0.SetLineStyle(2)
            line0.Draw("same")

            # legend=ROOT.TLegend(0.75,0.75,0.95,0.95)
            # legend.AddEntry(graph_ZLumi,"Z Lumi","pe")

            for graph in graphs_r:
                # legend.AddEntry(graph,graph.GetTitle(),"pe")
                graph.SetTitle("")
                graph.Draw("P same")

            # legend.AddEntry(graph_Lumi,refLumiSource+" Lumi","l")
            # legend.Draw("same")

            graph_rZLumi.Draw("P same")
            canvas.SaveAs(outDir+"/ZLumi"+str(fill)+"_"+xName+"_"+suffix+".pdf")

            canvas.SaveAs(subdir+"/ZLumi"+str(fill)+"_"+xName+"_"+suffix+".png")
            canvas.Close()
