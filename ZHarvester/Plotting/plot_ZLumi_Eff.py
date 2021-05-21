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
from ZUtils.python.utils import to_RootTime, getMCCorrection, cms, preliminary, workinprogress, text, custom_labels_y

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

parser.add_argument("-c", "--rates", required=True, type=str, help="csv file with z rates per measurement")
parser.add_argument("-s", "--saveDir",  default='./',  type=str, help="give output dir")
parser.add_argument("-f", "--fill",  default=0, type=int, help="specify a single fill to plot")
parser.add_argument("-m", "--mcCorrections",  default=None, type=str, help="specify MC correction file")
parser.add_argument("-y", "--year",  default="2017", type=str, help="specify year for MC comparison")
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

# >>>  load MC info and format
genInfoFile = 'res/GenInfo-V11_15-d20210427-t235949/infoMC_gen_{0}.json'.format(args.year)
recoInfoFile = 'res/RecoInfo-V11_15-d20210427-t235936/infoMC_reco_{0}.json'.format(args.year)

infoReco = pd.read_json(recoInfoFile,orient="index")
infoTrue = pd.read_json(genInfoFile,orient="index")
infoReco['bin']= infoReco.index
infoTrue['bin']= infoTrue.index

infoReco["True"] = False
infoTrue["True"] = True

info = pd.merge(infoTrue, infoReco, how='outer')
info = info.loc[info['bin'] != 'all']
info[['lo','hi']] = info['bin'].str.split("To", expand=True)
info['lo'] = info['lo'].apply(lambda x: int(x.replace("nPU","")))
info['hi'] = info['hi'].apply(lambda x: int(x.replace("Inf","75")))
info['center'] = info['lo'] + (info['hi'] - info['lo'])/2.

info.sort_values('center', inplace=True)

iReco = info.loc[info['True'] == False]
iTrue = info.loc[info['True'] == True]

iReco.reset_index(drop=True)
iTrue.reset_index(drop=True)

# pdb.set_trace()

########## Data Acquisition ##########

# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False) #, skiprows=[1,2,3,4,5])
if args.fill != 0:
    data = data.loc[data['fill'] == args.fill]

data = data.sort_values(['fill','tdate_begin','tdate_end'])

data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2
data['timeE'] = (data['tdate_end'] - data['tdate_begin'])/2

data['ZBBeff'] = data['ZBBeff'].apply(lambda x: unc.ufloat_fromstr(x))
data['ZBEeff'] = data['ZBEeff'].apply(lambda x: unc.ufloat_fromstr(x))
data['ZEEeff'] = data['ZEEeff'].apply(lambda x: unc.ufloat_fromstr(x))

if args.mcCorrections:
    print("Get MC corrections from file: "+args.mcCorrections)
    corr = getMCCorrection(args.mcCorrections)
    data['ZBBeff_mc'] = 1. / corr['effBB'](data['pileUp']) * data['ZBBeff']
    data['ZBEeff_mc'] = 1. / corr['effBE'](data['pileUp']) * data['ZBEeff']
    data['ZEEeff_mc'] = 1. / corr['effEE'](data['pileUp']) * data['ZEEeff']
else:
    data['ZBBeff_mc'] = data['ZBBeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))
    data['ZBEeff_mc'] = data['ZBEeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))
    data['ZEEeff_mc'] = data['ZEEeff_mc'].apply(lambda x: unc.ufloat_fromstr(x))

# seconds to hours
data['time'] = data['time'] / 3600.
data['timeE'] = data['timeE'] / 3600.
data['tdate_begin'] = data['tdate_begin'] / 3600.
data['tdate_end'] = data['tdate_end'] / 3600.

data = data.replace([np.inf, -np.inf], np.nan).dropna().dropna()


########## Plot ##########

for suffix, correction in (
    ("inclusive", 1 ),
    # ("BB",        'zLumiInstBB_mc', [("zLumiInstBB", 24, 2, "w/o MC"),] ),
    # ("BE",        'zLumiInstBE_mc', [("zLumiInstBE", 24, 2, "w/o MC"),] ),
    # ("EE",        'zLumiInstEE_mc', [("zLumiInstEE", 24, 2, "w/o MC"),] ),
    # ("All",       'All', [ ("zLumiInstBB_mc", 24, 2, "BB"), ("zLumiInstBE_mc", 25, 6, "BE"), ("zLumiInstEE_mc", 32, 4, "EE")]),
):
    ##### loop over Fills and produce fill specific plots
    for fill in data.drop_duplicates('fill')['fill'].values:
        dFill = data.loc[data['fill'] == fill]

        subdir = outDir+"/PlotsFill_"+str(fill)
        if not os.path.isdir(subdir):
            os.mkdir(subdir)

        starttime = dFill['time'].values[0]
        dFill_time = dFill['time'].values - starttime
        dFill_timeE = dFill['timeE'].values

        dFill_PU = dFill['pileUp'].values
        dFill_PUE = np.zeros(len(dFill_PU))

        if correction:
            yyEff = [dFill['ZBBeff_mc'].values,
                dFill['ZBEeff_mc'].values,
                dFill['ZEEeff_mc'].values]
        else:
            yyEff = [dFill['ZBBeff'].values,
                dFill['ZBEeff'].values,
                dFill['ZEEeff'].values]

        for yy, xxZ, xxZErr, xLabel, xName, legend_position in (
            (yyEff, dFill_time, dFill_timeE, "LHC Runtime [h]", "time", 1),
            (yyEff, dFill_PU,   dFill_PUE,   "<Number of pileup interactions>", "pileup", 2)
            ):

            canvas=ROOT.TCanvas("canvas","canvas",800,1000)
            pads = []
            graphs = []
            for i, l, padLo, padHi, margin_bottom, color, markerstyle in [
                (0, "BB",0.7, 1.0, 0.01, 2, 24),
                (1, "BE",0.4, 0.7, 0.01, 3, 25),
                (2, "EE",0.0, 0.4, 0.25, 4, 26)]:

                pad = ROOT.TPad("pad"+l, "pad "+l, 0., padLo, 1, padHi)
                canvas.SetTicks()
                pad.SetLeftMargin(margin_left)
                pad.SetRightMargin(margin_right)
                pad.SetTopMargin(margin_top)
                pad.SetBottomMargin(margin_bottom)
                pad.SetTickx()
                pad.SetTicky()
                pad.Draw()
                pad.cd()
                textsize1 = 28./(pad.GetWh()*pad.GetAbsHNDC())

                pads.append(pad)

                y = dFill['Z{0}eff_mc'.format(l)].values

                graph=ROOT.TGraphErrors(len(y), np.array(xxZ), np.array([iy.n for iy in y]), np.array(xxZErr), np.array([iy.s for iy in y]))
                # graph_Lumi.SetMarkerStyle(23)
                # graph_Lumi.SetMarkerColor(ROOT.kBlue+8)
                graph.SetMarkerStyle(markerstyle)
                graph.SetMarkerColor(color)
                graph.SetLineWidth(2)
                graph.SetMarkerSize(1.5)
                graph.SetTitle("")


                if xName == "pileup":
                    yTrue = iTrue['eff{0}'.format(l)].values
                    xTrue = iTrue['center'].values

                    gTrue=ROOT.TGraph(len(yTrue), xTrue, yTrue)
                    # gTrue.SetMarkerStyle(23)
                    # gTrue.SetMarkerColor(ROOT.kBlue+8)
                    gTrue.SetLineColor(color)
                    gTrue.SetLineWidth(2)
                    # gTrue.SetMarkerSize(1.5)
                    gTrue.SetTitle("")

                    graphs.append(gTrue)

                    yReco = iReco['eff{0}'.format(l)].values
                    xReco = iReco['center'].values

                    gReco=ROOT.TGraph(len(xReco), xReco, yReco)
                    # gTrue.SetMarkerStyle(23)
                    # gTrue.SetMarkerColor(ROOT.kBlue+8)
                    gReco.SetLineColor(color)
                    gReco.SetLineWidth(2)
                    gReco.SetLineStyle(2)
                    gReco.SetTitle("")

                    graphs.append(gReco)

                xmin = min((xxZ - xxZErr))
                xmax = max((xxZ + xxZErr))
                xWidth = xmax - xmin
                xmin = xmin - (xWidth * 0.05)
                xmax = xmax + (xWidth * 0.05)

                graph.GetYaxis().SetLabelSize(textsize1)

                if i == 1:
                    graph.GetYaxis().CenterTitle()
                    graph.GetYaxis().SetTitle(r"Z \rightarrow \mu\mu\ \mathrm{efficiency}")
                    graph.GetYaxis().SetTitleSize(textsize1)
                    graph.GetYaxis().SetTitleOffset(0.7)
                else:
                    graph.GetYaxis().SetTitle("")

                graph.GetYaxis().SetRangeUser(0.8,1.0)

                graph.GetYaxis().SetNdivisions(-204)
                # graph.GetYaxis().SetWmin(0.85)
                # graph.GetYaxis().SetWmax(0.95)

                custom_labels_y(graph, [" ", "0.85", "0.9", "0.95", " "])

                # graph_ZLumi.GetXaxis().SetTimeDisplay(1)
                graph.GetXaxis().SetTitleSize(textsize1)
                graph.GetXaxis().SetLabelSize(textsize1)
                graph.GetXaxis().SetTitle(xLabel)

                graph.GetXaxis().SetRangeUser(xmin,xmax)

                graphs.append(graph)

                graph.Draw("AP")

                if i == 2:
                    legend=ROOT.TLegend(0.13,0.26,0.975,0.44)
                else:
                    legend=ROOT.TLegend(0.13,0.06,0.975,0.24)
                # else:
                #     legend=ROOT.TLegend(0.84,0.06,0.975,0.24)

                legend.SetTextSize(textsize1)
                graphs.append(legend)

                legend.SetNColumns(3)

                legend.AddEntry(graph, l,"p")

                if i == 0 and xName == "pileup":

                    for labelname, linestyle in [("True MC",1), ("TnP MC",2)]:
                        entry = ROOT.TLegendEntry()
                        entry.SetLabel(labelname)
                        entry.SetOption("l")
                        entry.SetLineColor(ROOT.kGray)
                        entry.SetTextSize(textsize1)
                        entry.SetLineWidth(2)
                        entry.SetLineStyle(linestyle)
                        entry.SetTextAlign(12)
                        legend.GetListOfPrimitives().Add(entry)
                else:
                    for i in range(2):
                        legend.GetListOfPrimitives().Add(ROOT.TLegendEntry())


                graph.Draw("P same")

                if xName == "pileup":
                    gTrue.Draw("l same")
                    gReco.Draw("l same")
                legend.Draw("same")


                # cms(x=0.5, y=0.88, textsize=textsize1)
                if i == 0:
                    workinprogress(x=0.15, y=0.83, textsize=textsize1)
                    text("Fill "+str(fill), x=0.85, y=0.83, textsize=textsize1)

                graph.Draw("P same")
                canvas.cd()

            canvas.SaveAs(subdir+"/ZEfficiency"+str(fill)+"_"+xName+"_"+suffix+".eps")
            canvas.SaveAs(subdir+"/ZEfficiency"+str(fill)+"_"+xName+"_"+suffix+".png")
            canvas.Close()
