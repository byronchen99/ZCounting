import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
import shutil
import uncertainties as unc
import pdb

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cms", default="nothing", type=str, help="give the CMS csv as input")
parser.add_argument("-x", "--xsec", default="nothing", type=str, help="give the CMS csv as input")
parser.add_argument("-r", "--refLumi", default="nothing", type=str, help="give a ByLs.csv as input for reference Luminosity")
parser.add_argument("-s", "--saveDir", default='./', type=str, help="give output dir")
args = parser.parse_args()

if args.cms=="nothing":
	print "please provide cms input files"
	sys.exit()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Constants ##########

secPerLS=float(23.3)

MassMin_ = 56.
MassMax_ = 116.

#sigmaZ = 1870       # theory prediction from CMS PAS SMP-15-004
#acceptanceZ = 0.342367
#sigmaZfid = 610.1401700042975

#ZeffE=0.01          # systematic uncertainty of Z reconstruction efficiency
#sigmaZfidE = 0.03   # systematic uncertainty of fiducial Z cross section

### For plot labeling ###

currentYear = 2017

ptCut = 30
etaCut = 2.4
trigger='IsoMu27_v*'
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

data_xsec['zDelBB_mc'] = data_xsec['zDelBB_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data_xsec['zDelBE_mc'] = data_xsec['zDelBE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data_xsec['zDelEE_mc'] = data_xsec['zDelEE_mc'].apply(lambda x: unc.ufloat_fromstr(x))


xsecBB = sum(data_xsec['zDelBB_mc'])/sum(data_xsec['lumiRec'])
xsecBE = sum(data_xsec['zDelBE_mc'])/sum(data_xsec['lumiRec'])
xsecEE = sum(data_xsec['zDelEE_mc'])/sum(data_xsec['lumiRec'])

xsec = xsecBB + xsecBE + xsecEE

# --- z luminosity
data = pd.read_csv(str(args.cms), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(['fill','tdate_begin','tdate_end'])


data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2
data['timeE'] = (data['tdate_end'] - data['tdate_begin'])/2

data['zDelBB_mc'] = data['zDelBB_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDelBE_mc'] = data['zDelBE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDelEE_mc'] = data['zDelEE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDel_mc'] = data['zDelBB_mc'] + data['zDelBE_mc'] + data['zDelEE_mc']
data['zRate_mc'] = data['zDel_mc'] / data['timewindow']
data['zLumiInst_mc'] = data['zDel_mc'] / (data['timewindow'] * xsec)

data['zXSec_mc'] = data['zDel_mc'] / data['lumiDel']

zcountlist = data.groupby('fill')['zDel_mc'].apply(list)

#sort out points with low statistics
#data = data[data['z_relstat'] < 0.05]

data = data.replace([np.inf, -np.inf], np.nan).dropna().dropna()

zcountl = []
timel = []

suffix=""

if "Central" in str(args.cms):
	suffix="Barrel"
else:
	suffix="Inclusive"

metaFills= []
metaXsecCMS= []

slope = dict()
chi2 = dict()

zcountsAccu=0
metazcountsAccu=[]
metazcountsoverlumi=[]

########## Plot ##########

### chi2 values of each category
# for c in ('HLTeffB_chi2pass', 'HLTeffB_chi2fail', 'HLTeffE_chi2pass', 'HLTeffE_chi2fail',
#           'SeleffB_chi2pass', 'SeleffB_chi2fail', 'SeleffE_chi2pass', 'SeleffE_chi2fail',
#           # 'GloeffB_chi2pass', 'GloeffB_chi2fail', 'GloeffE_chi2pass', 'GloeffE_chi2fail',
#           'TrkeffB_chi2pass', 'TrkeffB_chi2fail', 'TrkeffE_chi2pass', 'TrkeffE_chi2fail',
#           'StaeffB_chi2pass', 'StaeffB_chi2fail', 'StaeffE_chi2pass', 'StaeffE_chi2fail',
#          ):
#
#     graph_chi2 = ROOT.TGraph(len(data), data['time'].values, data[c].values)
#     graph_chi2.SetName("graph_chi2")
#     graph_chi2.SetMarkerStyle(23)
#     graph_chi2.SetMarkerColor(ROOT.kAzure-4)
#     graph_chi2.SetMarkerSize(1.5)
#     graph_chi2.SetTitle(c)
#
#     graph_chi2.GetYaxis().SetTitle("#chi^{2}/ndf")
#     graph_chi2.GetXaxis().SetTitle("Time")
#     graph_chi2.GetXaxis().SetTimeDisplay(1)
#     graph_chi2.GetXaxis().SetTimeOffset(0,"gmt")
#     graph_chi2.GetXaxis().SetTitleSize(0.06)
#     graph_chi2.GetYaxis().SetTitleSize(0.06)
#     graph_chi2.GetXaxis().SetTitleOffset(0.72)
#     graph_chi2.GetYaxis().SetTitleOffset(1.1)
#     graph_chi2.GetXaxis().SetLabelSize(0.05)
#     graph_chi2.GetYaxis().SetLabelSize(0.05)
#     #graph_chi2.GetYaxis().SetRangeUser(-0.01,0.01)
#     c3=ROOT.TCanvas("c3_"+c,"c3 "+c,1000,600)
#     c3.SetGrid()
#
#     # mean, where outlier with sigma > 1 are rejected
#     avg_chi2 = np.mean(data[c][abs(data[c] - np.mean(data[c])) < np.std(data[c])])
#
#     graph_chi2.Draw("AP")
#
#     legend=ROOT.TLegend(0.2,0.6,0.4,0.7)
#     legend.AddEntry(graph_chi2,"Measurement","p")
#     legend.Draw("same")
#
#     #text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
#     #text.SetNDC()
#     #text.Draw()
#     text2=ROOT.TLatex(0.2,0.23,"avg chi2: "+str(avg_chi2))
#     text2.SetNDC()
#     text2.Draw()
#
#     c3.SaveAs(outDir+"/"+c+".png")
#     c3.Close()
#
#     th1_chi2 = ROOT.TH1F("h1_"+c, "h1 title "+c, 20, 0., 2.0*avg_chi2)
#     th1_chi2.GetXaxis().SetTitle("#chi^{2}/ndf")
#     th1_chi2.GetYaxis().SetTitle("fits")
#     th1_chi2.GetXaxis().SetTitleSize(0.06)
#     th1_chi2.GetYaxis().SetTitleSize(0.06)
#     th1_chi2.GetXaxis().SetTitleOffset(0.72)
#     th1_chi2.GetYaxis().SetTitleOffset(1.1)
#     th1_chi2.GetXaxis().SetLabelSize(0.05)
#     th1_chi2.GetYaxis().SetLabelSize(0.05)
#     th1_chi2.SetTitle(c)
#     th1_chi2.SetLineColor(1)
#     th1_chi2.SetLineWidth(2)
#     th1_chi2.SetFillColor(0)
#
#     for v in data[c].values:
#         th1_chi2.Fill(v)
#
#     #f_chi2 = ROOT.TF1("fchi2","{0}*TMath::Prob(x,1)".format(th1_chi2.Integral()*(2.0*avg_chi2)/20.),0.0,2.0*avg_chi2)
#     #f_chi2.SetLineWidth(2)
#     #f_chi2.SetLineColor(2)
#
#     #legend=ROOT.TLegend(0.8,0.8,0.9,0.9)
#     #legend.AddEntry(f_chi2,"chi2(x,1)","l")
#     #legend.AddEntry(th1_chi2,"fits","f")
#
#     c4=ROOT.TCanvas("c4_"+c,"c4_"+c,1000,600)
#
#     th1_chi2.Draw('HIST')
#     #f_chi2.Draw("same")
#     #legend.Draw("same")
#
#
#     c4.SaveAs(outDir+"/hist_"+c+".png")
#     c4.Close()



##### loop over Fills and produce fill specific plots
for fill in data.drop_duplicates('fill')['fill'].values:
    dFill = data.loc[data['fill'] == fill]

    zcountl.append(sum(zcountlist[fill]))
    zcountsAccu=zcountsAccu+sum(zcountlist[fill])
    metazcountsAccu.append(zcountsAccu)
    timel.append(dFill.iloc[-1]['time'])

    startTime=dFill.iloc[0]['time']
    endTime = dFill.iloc[-1]['time']

    shutil.rmtree(outDir+"/PlotsFill_"+str(fill), ignore_errors=True)
    os.makedirs(outDir+"/PlotsFill_"+str(fill))

    ### Efficiency vs Pileup plots ###

    for eff, name in (
        # ('ZBBeff_mc','corrected Z-BB-Reconstruction efficitency'),
        # ('ZBEeff_mc','corrected Z-BE-Reconstruction efficitency'),
        # ('ZEEeff_mc','corrected Z-EE-Reconstruction efficitency'),
        ('ZBBeff'  ,'Z-BB-Reconstruction efficitency'),
        ('ZBEeff'  ,'Z-BE-Reconstruction efficitency'),
        ('ZEEeff'  ,'Z-EE-Reconstruction efficitency'),
        ('HLTeffB' ,'Muon HLT-B efficiency'),
        ('HLTeffE' ,'Muon HLT-E efficiency'),
        ('SeleffB' ,'Muon Sel-B efficiency'),
        ('SeleffE' ,'Muon Sel-E efficiency'),
        # ('GloeffB' ,'Muon Glo-B efficiency'),
        # ('GloeffE' ,'Muon Glo-E efficiency'),
        ('TrkeffB' ,'Muon Trk-B efficiency'),
        ('TrkeffE' ,'Muon Trk-E efficiency'),
        ('StaeffB' ,'Muon Sta-B efficiency'),
        ('StaeffE' ,'Muon Sta-E efficiency'),
    ):

        yy = dFill[eff].apply(lambda x: unc.ufloat_fromstr(x).nominal_value).values
        yyErr = dFill[eff].apply(lambda x: unc.ufloat_fromstr(x).std_dev).values

        graph_Zeff = ROOT.TGraphErrors(len(dFill), dFill['pileUp'].values, yy, yyErr)
        graph_Zeff.SetName("graph_Zeff")
        graph_Zeff.SetMarkerStyle(22)
        graph_Zeff.SetMarkerColor(ROOT.kOrange+8)
        graph_Zeff.SetFillStyle(0)
        graph_Zeff.SetMarkerSize(1.5)
        graph_Zeff.SetTitle(name+", Fill "+str(fill))
        graph_Zeff.GetYaxis().SetTitle(eff)
        graph_Zeff.GetYaxis().SetTitleSize(0.07)
        graph_Zeff.GetYaxis().SetTitleOffset(1.1)
        graph_Zeff.GetXaxis().SetTitle("avg. PU")
        graph_Zeff.GetXaxis().SetTitleSize(0.06)
        graph_Zeff.GetXaxis().SetTitleOffset(0.75)
        graph_Zeff.GetXaxis().SetLabelSize(0.05)
        graph_Zeff.GetYaxis().SetLabelSize(0.05)

        fr = graph_Zeff.Fit("pol1","S")

        if eff in slope.keys():
            slope[eff].append(fr.Parameter(1))
        else:
            slope[eff] = [fr.Parameter(1)]

        c1=ROOT.TCanvas("c1","c1",1000,600)
        c1.cd(1)
        graph_Zeff.Draw("AP")

        legend=ROOT.TLegend(0.65,0.65,0.9,0.9)
        legend.AddEntry(graph_Zeff,"CMS","pe")

        #text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        #text1.SetNDC()
        #text1.Draw()

        text_fr =  ROOT.TText(0.2,0.30,"p0 = "+str(round(fr.Parameter(0),5))+" \pm "+str(round(fr.ParError(0),5)) )
        text_fr1 = ROOT.TText(0.2,0.23,"p1 = "+str(round(fr.Parameter(1),5))+" \pm "+str(round(fr.ParError(1),5)))
        text_fr2 = ROOT.TText(0.2,0.16,"Chi2: "+str(round(fr.Chi2(),5))+" Ndf: "+str(fr.Ndf()))
        text_fr.SetNDC()
        text_fr.Draw()
        text_fr1.SetNDC()
        text_fr1.Draw()
        text_fr2.SetNDC()
        text_fr2.Draw()

        c1.SaveAs(outDir+"/PlotsFill_"+str(fill)+"/"+str(eff)+"_PU"+str(fill)+".png")
        c1.Close()


    ### Z Rates ###

    yy = dFill['zRate_mc'].apply(unc.nominal_value).values
    yyErr = dFill['zRate_mc'].apply(unc.std_dev).values

    graph_zrates=ROOT.TGraphErrors(len(dFill),dFill['time'].values, yy, dFill['timeE'].values, yyErr )
    graph_zrates.SetName("graph_zrates")
    graph_zrates.SetMarkerStyle(22)
    graph_zrates.SetMarkerColor(ROOT.kOrange+8)
    graph_zrates.SetFillStyle(0)
    graph_zrates.SetMarkerSize(1.5)
    graph_zrates.GetXaxis().SetTimeDisplay(1)
    graph_zrates.SetTitle(suffix+" Z-Rates, Fill "+str(fill))
    graph_zrates.GetYaxis().SetTitle("Z-Rate [Hz]")
    graph_zrates.GetYaxis().SetTitleSize(0.07)
    graph_zrates.GetYaxis().SetTitleOffset(0.7)
    graph_zrates.GetXaxis().SetTitle("Time")
    graph_zrates.GetXaxis().SetTitleSize(0.06)
    graph_zrates.GetXaxis().SetTitleOffset(0.75)
    graph_zrates.GetXaxis().SetLabelSize(0.05)
    graph_zrates.GetYaxis().SetLabelSize(0.05)
    graph_zrates.GetXaxis().SetRangeUser(startTime,endTime)

    c1=ROOT.TCanvas("c1","c1",1000,600)
    c1.cd(1)
    graph_zrates.Draw("AP")

    legend=ROOT.TLegend(0.65,0.65,0.9,0.9)
    legend.AddEntry(graph_zrates,"CMS","pe")

    #text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text1.SetNDC()
    #text1.Draw()
    c1.SaveAs(outDir+"/PlotsFill_"+str(fill)+"/zrates"+str(fill)+suffix+".png")
    c1.Close()

    ### Z Luminosity vs reference Luminosity ###

    dFill_ref = data_ref.loc[data_ref['fill'] == fill]
    graph_Lumi=ROOT.TGraph(len(dFill_ref),dFill_ref['time'].values,  dFill_ref['dLRec(/pb)'].values)
    graph_Lumi.SetName("graph_Lumi")
    graph_Lumi.SetMarkerStyle(23)
    graph_Lumi.SetMarkerColor(ROOT.kBlue+8)
    graph_Lumi.SetFillStyle(0)
    graph_Lumi.SetMarkerSize(1.5)

    yy = dFill['zLumiInst_mc'].apply(unc.nominal_value).values
    yyErr = dFill['zLumiInst_mc'].apply(unc.std_dev).values

    graph_ZLumi=ROOT.TGraphErrors(len(dFill),dFill['time'].values, yy, dFill['timeE'].values, yyErr )
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
    c2.SaveAs(outDir+"/PlotsFill_"+str(fill)+"/ZLumi"+str(fill)+suffix+".png")
    c2.Close()

    ### Cross sections ###

    metaXsecCMS.append(sum(dFill['zXSec_mc'].values)/len(dFill))
    metaFills.append(float(fill))

    yy = dFill['zXSec_mc'].apply(unc.nominal_value).values
    yyErr = dFill['zXSec_mc'].apply(unc.std_dev).values

    graph_xsec=ROOT.TGraphErrors(len(dFill), dFill['time'].values, yy, yyErr)
    graph_xsec.SetName("graph_xsec")
    graph_xsec.SetMarkerStyle(22)
    graph_xsec.SetMarkerColor(ROOT.kOrange+8)
    graph_xsec.SetFillStyle(0)
    graph_xsec.SetMarkerSize(1.5)
    graph_xsec.SetFillStyle(0)
    graph_xsec.SetTitle(" cross section, Fill "+str(fill))
    graph_xsec.GetXaxis().SetTimeDisplay(1)
    graph_xsec.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
    graph_xsec.GetYaxis().SetTitleSize(0.06)
    graph_xsec.GetYaxis().SetTitleOffset(0.80)
    graph_xsec.GetXaxis().SetTitle("Time")
    graph_xsec.GetXaxis().SetTitleSize(0.06)
    graph_xsec.GetXaxis().SetTitleOffset(0.72)
    graph_xsec.GetXaxis().SetLabelSize(0.05)
    graph_xsec.GetYaxis().SetLabelSize(0.05)

    c4=ROOT.TCanvas("c4","c4",1000,600)
    c4.SetGrid()
    graph_xsec.Draw("AP")

    legend=ROOT.TLegend(0.75,0.75,0.9,0.9)
    legend.AddEntry(graph_xsec,"CMS","p")

    #text=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text.SetNDC()
    #text.Draw()
    text2=ROOT.TLatex(0.6,0.23,"#splitline{"+str(MassMin_)+" GeV<M(#mu#mu) < "+str(MassMax_)+" GeV}{p_{T}(#mu)>"+str(ptCut)+" GeV, |#eta(#mu)|<"+str(etaCut)+"}")
    text2.SetNDC()
    text2.SetTextSize(0.04)
    text2.Draw()
    c4.SaveAs(outDir+"/PlotsFill_"+str(fill)+"/ZStability"+str(fill)+suffix+".png")

    c4.Close()


### Efficiency-vs-Pileup slopes vs fill
for key, values in slope.items():

    i_slope = np.array(values)

    graph_slopesEffPu = ROOT.TGraph(len(metaFills), np.array(metaFills), i_slope)
    graph_slopesEffPu.SetName("graph_slopesEffPu")
    graph_slopesEffPu.SetMarkerStyle(22)
    graph_slopesEffPu.SetMarkerColor(ROOT.kOrange+8)
    graph_slopesEffPu.SetMarkerSize(2.5)
    graph_slopesEffPu.SetTitle(key+" efficiency-vs-pileup slope")

    graph_slopesEffPu.GetXaxis().SetTitle("Fill")
    graph_slopesEffPu.GetYaxis().SetTitle("d"+key+"/dPU")
    graph_slopesEffPu.GetXaxis().SetTitleSize(0.06)
    graph_slopesEffPu.GetYaxis().SetTitleSize(0.06)
    graph_slopesEffPu.GetXaxis().SetTitleOffset(0.72)
    graph_slopesEffPu.GetYaxis().SetTitleOffset(1.1)
    graph_slopesEffPu.GetXaxis().SetLabelSize(0.05)
    graph_slopesEffPu.GetYaxis().SetLabelSize(0.05)
    graph_slopesEffPu.GetYaxis().SetRangeUser(-0.01,0.01)
    c3=ROOT.TCanvas("c3","c3",1000,600)
    c3.SetGrid()

    # mean, where outlier with sigma > 1 are rejected
    avg_slope = np.mean(i_slope[abs(i_slope - np.mean(i_slope)) < np.std(i_slope)])

    graph_slopesEffPu.Draw("AP")

    #text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text.SetNDC()
    #text.Draw()
    text2=ROOT.TLatex(0.2,0.23,"avg slope: "+str(avg_slope))
    text2.SetNDC()
    text2.Draw()

    c3.SaveAs(outDir+"/slopes"+key+"vsPU.png")
    c3.Close()


### fiducial cross section vs fill

ROOT.gROOT.SetBatch(True)
metaXsecCMS2= []
for n in range(0,len(metaXsecCMS)):
	metaXsecCMS2.append(metaXsecCMS[n]/(sum(metaXsecCMS)/len(metaXsecCMS)))


graph_metacmsXsec=ROOT.TGraph(len(metaFills), np.array(metaFills), np.array(metaXsecCMS))
graph_metacmsXsec.SetName("graph_metaXsecCms")
graph_metacmsXsec.SetMarkerStyle(22)
graph_metacmsXsec.SetMarkerColor(ROOT.kOrange+8)
graph_metacmsXsec.SetMarkerSize(2.5)
graph_metacmsXsec.SetTitle("Cross Section Summary, "+suffix+" Z-Rates")

multMetaGraphXsec=ROOT.TMultiGraph("multMetaGraphXsec", suffix+" Z-Rates")
multMetaGraphXsec.SetName("multMetaGraphXsec")
graph_metacmsXsec.GetXaxis().SetTitle("Fill")
graph_metacmsXsec.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
graph_metacmsXsec.GetXaxis().SetTitleSize(0.06)
graph_metacmsXsec.GetYaxis().SetTitleSize(0.06)
graph_metacmsXsec.GetXaxis().SetTitleOffset(0.72)
graph_metacmsXsec.GetYaxis().SetTitleOffset(0.8)
graph_metacmsXsec.GetXaxis().SetLabelSize(0.05)
graph_metacmsXsec.GetYaxis().SetLabelSize(0.05)

multMetaGraphXsec.Add(graph_metacmsXsec)

c3=ROOT.TCanvas("c3","c3",1000,600)
c3.SetGrid()

graph_metacmsXsec.Draw("AP")


#text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
#text.SetNDC()
#text.Draw()
text2=ROOT.TLatex(0.6,0.23,"#splitline{"+str(MassMin_)+" GeV<M(#mu#mu) < "+str(MassMax_)+" GeV}{p_{T}(#mu)>"+str(ptCut)+" GeV, |#eta(#mu)|<"+str(etaCut)+"}")
text2.SetNDC()
text2.SetTextSize(0.04)
text2.Draw()
c3.SaveAs(outDir+"/summaryZStability"+suffix+".png")
c3.Close()

### corrected Z Counts vs fill
yy = np.array([y.n for y in zcountl])
yyErr = np.array([y.s for y in zcountl])

graph_zcount=ROOT.TGraphErrors(len(metaFills), np.array(metaFills), yy, yyErr)
graph_zcount.SetName("graph_zcount")
graph_zcount.SetMarkerStyle(22)
graph_zcount.SetMarkerColor(ROOT.kOrange+8)
graph_zcount.SetMarkerSize(2.5)
graph_zcount.SetTitle("Z Counts Per Fill")
graph_zcount.GetXaxis().SetTitle("Fill")
graph_zcount.GetYaxis().SetTitle("Z Count")
graph_zcount.GetXaxis().SetTitleSize(0.06)
graph_zcount.GetYaxis().SetTitleSize(0.06)
graph_zcount.GetXaxis().SetTitleOffset(0.72)
graph_zcount.GetYaxis().SetTitleOffset(0.8)
graph_zcount.GetXaxis().SetLabelSize(0.05)
graph_zcount.GetYaxis().SetLabelSize(0.05)
graph_zcount.GetXaxis().SetLabelSize(0.05)
graph_zcount.GetYaxis().SetLabelSize(0.05)


c5=ROOT.TCanvas("c5","c5",1000,600)
c5.SetGrid()
graph_zcount.Draw("AP")

#text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
#text.SetNDC()
#text.Draw()
text2=ROOT.TLatex(0.6,0.23,"#splitline{"+str(MassMin_)+" GeV<M(#mu#mu) < "+str(MassMax_)+" GeV}{p_{T}(#mu)>"+str(ptCut)+" GeV, |#eta(#mu)|<"+str(etaCut)+"}")
text2.SetNDC()
text2.SetTextSize(0.04)
text2.Draw()
text3=ROOT.TLatex(0.6,0.33,"#color[4]{"+trigger+"}")
text3.SetNDC()
text3.SetTextSize(0.04)
text3.Draw()
c5.SaveAs(outDir+"/ZCountPerFill"+suffix+".png")
c5.Close()



### Accumalated corrected Z counts
yy = np.array([y.n for y in metazcountsAccu])
yyErr = np.array([y.s for y in metazcountsAccu])

graph_zcountA=ROOT.TGraphErrors(len(metaFills), np.array(timel), yy, yyErr)
graph_zcountA.SetName("graph_zcountAccu")
graph_zcountA.SetMarkerStyle(22)
graph_zcountA.SetMarkerColor(ROOT.kOrange+8)
graph_zcountA.SetMarkerSize(2.5)
graph_zcountA.SetTitle("Accumulated Z Bosons over Time")
graph_zcountA.GetXaxis().SetTitle("Time")
graph_zcountA.GetXaxis().SetTimeDisplay(1)
graph_zcountA.GetXaxis().SetTimeOffset(0,"gmt")
graph_zcountA.GetYaxis().SetTitle("Z Count")
graph_zcountA.GetXaxis().SetTitleSize(0.06)
graph_zcountA.GetYaxis().SetTitleSize(0.06)
graph_zcountA.GetXaxis().SetTitleOffset(0.72)
graph_zcountA.GetYaxis().SetTitleOffset(0.8)
graph_zcountA.GetXaxis().SetLabelSize(0.05)
graph_zcountA.GetYaxis().SetLabelSize(0.05)
graph_zcountA.GetXaxis().SetLabelSize(0.05)
graph_zcountA.GetYaxis().SetLabelSize(0.05)

c6=ROOT.TCanvas("c6","c6",1000,600)
c6.SetGrid()
graph_zcountA.Draw("AP")
#text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
#text.SetNDC()
#text.Draw()
text2=ROOT.TLatex(0.6,0.23,"#splitline{"+str(MassMin_)+" GeV<M(#mu#mu) < "+str(MassMax_)+" GeV}{p_{T}(#mu)>"+str(ptCut)+" GeV, |#eta(#mu)|<"+str(etaCut)+"}")
text2.SetNDC()
text2.SetTextSize(0.04)
text2.Draw()
text3=ROOT.TLatex(0.6,0.33,"#color[4]{"+trigger+"}")
text3.SetNDC()
text3.SetTextSize(0.04)
text3.Draw()
c6.SaveAs(outDir+"/ZCountAccumulated"+suffix+".png")
c6.Close()
