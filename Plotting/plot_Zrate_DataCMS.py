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
parser.add_argument("-r", "--refLumi", default="nothing", type=str, help="give a ByLs.csv as input for reference Luminosity")
parser.add_argument("-s", "--saveDir", default='./', type=str, help="give output dir")
args = parser.parse_args()

if args.cms=="nothing":
	print "please provide cms input files"
	sys.exit()

print args.cms

########## Constants ##########

secPerLS=float(23.3)

sigmaZ = 1870       # theory prediction from CMS PAS SMP-15-004  
acceptanceZ = 0.342367
sigmaZfid = sigmaZ * acceptanceZ

ZeffE=0.03          # systematic uncertainty of Z reconstruction efficiency
sigmaZfidE = 0.05   # systematic uncertainty of fiducial Z cross section

### For plot labeling ###

ptCut = 27
etaCut = 2.4
trigger='IsoMu24_v*'
refLumiSource='hfoc18v5'

########## Data Acquisition ##########

data = pandas.read_csv(str(args.cms), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(['fill','beginTime','endTime'])
zcountlist = data.groupby('fill')['delZCount'].apply(list)
delLumilist = data.groupby('fill')['delLumi'].apply(list)
timelist = data.groupby('fill')['endTime'].apply(list)

# convert reference data from ByLs.csv file
lumiFile=open(str(args.refLumi))
lumiLines=lumiFile.readlines()
data_ref = pandas.read_csv(str(args.refLumi), sep=',',low_memory=False, skiprows=lambda x: lumiLines[x].startswith('#') and not lumiLines[x].startswith('#run'))
if 'delivered(/ub)' in data_ref.columns.tolist():      #convert to /pb
    data_ref['delivered(/ub)'] = data_ref['delivered(/ub)'].apply(lambda x:x / 1000000.)
    data_ref['recorded(/ub)'] = data_ref['recorded(/ub)'].apply(lambda x:x / 1000000.)
    data_ref = data_ref.rename(index=str, columns={'delivered(/ub)':'delivered(/pb)', 'recorded(/ub)':'recorded(/pb)' })

data_ref['fill'] = pandas.to_numeric(data_ref['#run:fill'].str.split(':',expand=True)[1])
data_ref['run'] = pandas.to_numeric(data_ref['#run:fill'].str.split(':',expand=True)[0])
data_ref['ls'] = pandas.to_numeric(data_ref['ls'].str.split(':',expand=True)[0])
data_ref = data_ref.drop(['#run:fill','hltpath','source'],axis=1)
data_ref = data_ref.sort_values(['fill','run','ls','delivered(/pb)','recorded(/pb)'])
data_ref = data_ref.drop_duplicates(['fill','run','ls'])

data_ref[['date','time']] = data_ref['time'].str.split(" ",expand=True)
data_ref[['month','day','year']] = (data_ref['date'].str.split("/",expand=True))
data_ref[['hour','min','sec']] = (data_ref['time'].str.split(":",expand=True))
data_ref['tdate'] = "20"+data_ref['year']+"-"+data_ref['month']+"-"+data_ref['day']+" "+data_ref['hour']+":"+data_ref['min']+":"+data_ref['sec']
data_ref['tdate'] = data_ref['tdate'].apply(ROOT.TDatime)
data_ref['tdate'] = data_ref['tdate'].apply(lambda x: x.Convert()).astype(float)
data_ref['tdateE'] = data_ref['tdate'] - data_ref['tdate']
data_ref['dLdel(/pb)'] = data_ref['delivered(/pb)']/secPerLS 	#delivered instantanious luminosity
data_ref = data_ref[['fill','dLdel(/pb)','tdate']]		#Keep only what you need 


# remove rows with invalid Z rates and low statistics
data = data[np.isfinite(data['ZRate'])]
data = data[data['delZCount'] > 1000.]

# convert time into root TDatime format
data[['date','time']] = data['beginTime'].str.split(" ",expand=True)
data[['month','day','year']] = (data['date'].str.split("/",expand=True))
data[['hour','min','sec']] = (data['time'].str.split(":",expand=True))
data['tdate'] = "20"+data['year']+"-"+data['month']+"-"+data['day']+" "+data['hour']+":"+data['min']+":"+data['sec']
data['tdate'] = data['tdate'].apply(ROOT.TDatime)
data['tdate'] = data['tdate'].apply(lambda x: x.Convert()).astype(float)
data['tdateE'] = data['tdate'] - data['tdate']
data = data.drop(['date','time','month','day','year','hour','min','sec'],axis=1)

data['ZRateE'] = data['ZRate'] * np.sqrt( 1./data['delZCount'] + ZeffE**2 )

data['Xsec'] = data['delZCount']/data['delLumi']
data['XsecEy'] = data['ZRate']/data['instDelLumi']*0.02

data['ZLumi'] = data['ZRate']/sigmaZfid
data['ZLumiE'] = data['ZLumi'] * np.sqrt((data['ZRateE']/data['ZRate'])**2 + sigmaZfidE**2 )

data['instDelLumiE'] = data['instDelLumi']*0.02

zcountl=array('d')
timel=array('d')

suffix=""

if "Central" in str(args.cms):
	suffix="Barrel"
else:
	suffix="Inclusive"

metaFills=array('d')
metaXsecCMS=array('d')

zcountsAccu=0
metazcountsAccu=array('d')
metazcountsoverlumi=array('d')

########## Plot ##########

for fill in data.drop_duplicates('fill')['fill'].values:
	dFill = data.loc[data['fill'] == fill]
        
        zcountl.append(sum(zcountlist[fill]))
	zcountsAccu=zcountsAccu+sum(zcountlist[fill])
	metazcountsAccu.append(zcountsAccu)
	timel.append(dFill.iloc[-1]['tdate'])
  
	
	startTime=dFill.iloc[0]['tdate']
	endTime = dFill.iloc[-1]['tdate']
	
        shutil.rmtree(args.saveDir+"PlotsFill_"+str(fill), ignore_errors=True)
	os.makedirs(args.saveDir+"PlotsFill_"+str(fill))
	
        ### Z Rates ###
        
        graph_cms=ROOT.TGraphErrors(len(dFill),dFill['tdate'].values,dFill['ZRate'].values, dFill['tdateE'].values, dFill['ZRateE'].values )
        graph_cms.SetName("graph_cms")
        graph_cms.SetMarkerStyle(22)
        graph_cms.SetMarkerColor(ROOT.kOrange+8)
        graph_cms.SetFillStyle(0)
        graph_cms.SetMarkerSize(1.5)
	graph_cms.GetXaxis().SetTimeDisplay(1)
	graph_cms.SetTitle(suffix+" Z-Rates, Fill "+str(fill))
	graph_cms.GetYaxis().SetTitle("Z-Rate [Hz]")
	graph_cms.GetYaxis().SetTitleSize(0.07)
	graph_cms.GetYaxis().SetTitleOffset(0.7)
        graph_cms.GetXaxis().SetTitle("Time")
	graph_cms.GetXaxis().SetTitleSize(0.06)
	graph_cms.GetXaxis().SetTitleOffset(0.75)
	graph_cms.GetXaxis().SetLabelSize(0.05)
	graph_cms.GetYaxis().SetLabelSize(0.05)
	graph_cms.GetXaxis().SetRangeUser(startTime,endTime)
		
        c1=ROOT.TCanvas("c1","c1",1000,600)
	c1.cd(1)
	graph_cms.Draw("AP")

	legend=ROOT.TLegend(0.65,0.65,0.9,0.9)
	legend.AddEntry(graph_cms,"CMS","p")

	text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
	text1.SetNDC()
	text1.Draw()
	c1.SaveAs(args.saveDir+"PlotsFill_"+str(fill)+"/zrates"+str(fill)+suffix+".png")
	c1.Close()

        ### Z Luminosity vs reference Luminosity ###

        dFill_ref = data_ref.loc[data_ref['fill'] == fill]
        graph_Lumi=ROOT.TGraph(len(dFill_ref),dFill_ref['tdate'].values,  dFill_ref['dLdel(/pb)'].values)
        graph_Lumi.SetName("graph_Lumi")
        graph_Lumi.SetMarkerStyle(23)
        graph_Lumi.SetMarkerColor(ROOT.kBlue+8)
        graph_Lumi.SetFillStyle(0)
        graph_Lumi.SetMarkerSize(1.5)

        graph_ZLumi=ROOT.TGraphErrors(len(dFill),dFill['tdate'].values,dFill['ZLumi'].values, dFill['tdateE'].values, dFill['ZLumiE'].values )
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
        graph_ZLumi.GetXaxis().SetRangeUser(startTime,endTime)

        graph_ZLumi.SetTitle(suffix+" inst. Luminosity, Fill "+str(fill))

        c2=ROOT.TCanvas("c2","c2",1000,600)
        c2.cd(1)
        graph_ZLumi.Draw("AP")
        graph_Lumi.Draw("l same")

        legend=ROOT.TLegend(0.65,0.65,0.9,0.8)
        legend.AddEntry(graph_ZLumi,"Z Lumi","p")
        legend.AddEntry(graph_Lumi,refLumiSource+" Lumi","l")
        legend.Draw("same")

        text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        text1.SetNDC()
        text1.Draw()
        c2.SaveAs(args.saveDir+"PlotsFill_"+str(fill)+"/ZLumi"+str(fill)+suffix+".png")
        c2.Close()

        ### Cross sections ###

	metaXsecCMS.append(sum(dFill['Xsec'].values)/len(dFill))	
	metaFills.append(float(fill))	

	graph_cmsXsec=ROOT.TGraph(len(dFill),dFill['tdate'].values,dFill['Xsec'].values)
        graph_cmsXsec.SetName("graph_cmsXsec")
        graph_cmsXsec.SetMarkerStyle(22)
        graph_cmsXsec.SetMarkerColor(ROOT.kOrange+8)
        graph_cmsXsec.SetFillStyle(0)
        graph_cmsXsec.SetMarkerSize(1.5)
	graph_cmsXsec.SetFillStyle(0)
	graph_cmsXsec.SetTitle(" cross section, Fill "+str(fill))
        graph_cmsXsec.GetXaxis().SetTimeDisplay(1)
	graph_cmsXsec.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
	graph_cmsXsec.GetYaxis().SetTitleSize(0.06)
	graph_cmsXsec.GetYaxis().SetTitleOffset(0.80)
	graph_cmsXsec.GetXaxis().SetTitle("Time")
	graph_cmsXsec.GetXaxis().SetTitleSize(0.06)
	graph_cmsXsec.GetXaxis().SetTitleOffset(0.72)
	graph_cmsXsec.GetXaxis().SetLabelSize(0.05)
	graph_cmsXsec.GetYaxis().SetLabelSize(0.05)	

	c4=ROOT.TCanvas("c4","c4",1000,600)
	c4.SetGrid()	
	graph_cmsXsec.Draw("AP")		
	
	legend=ROOT.TLegend(0.75,0.75,0.9,0.9)
	legend.AddEntry(graph_cmsXsec,"CMS","p")
	
        text=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
	text.SetNDC()
	text.Draw()
	text2=ROOT.TLatex(0.6,0.23,"#splitline{66 GeV<M(#mu#mu) < 116 GeV}{p_{T}(#mu)>"+str(ptCut)+" GeV, |#eta(#mu)|<"+str(etaCut)+"}")
	text2.SetNDC()
	text2.SetTextSize(0.04)
	text2.Draw()
	c4.SaveAs(args.saveDir+"PlotsFill_"+str(fill)+"/ZStability"+str(fill)+suffix+".png")
		
	c4.Close()
	

### fiducial cross section per fill
	
ROOT.gROOT.SetBatch(True)
metaXsecCMS2=array('d')
for n in range(0,len(metaXsecCMS)):
	metaXsecCMS2.append(metaXsecCMS[n]/(sum(metaXsecCMS)/len(metaXsecCMS)))	


graph_metacmsXsec=ROOT.TGraph(len(metaFills),metaFills,metaXsecCMS)
graph_metacmsXsec.SetName("graph_metaXsecCms")
graph_metacmsXsec.SetMarkerStyle(22)
graph_metacmsXsec.SetMarkerColor(ROOT.kOrange+8)
graph_metacmsXsec.SetMarkerSize(2.5)
graph_metacmsXsec.SetTitle("Cross Section Summary, "+suffix+" Z-Rates")

multMetaGraphXsec=ROOT.TMultiGraph("multMetaGraphXsec",suffix+" Z-Rates")
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


text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
text.SetNDC()
text.Draw()
text2=ROOT.TLatex(0.6,0.23,"#splitline{66 GeV<M(#mu#mu) < 116 GeV}{p_{T}(#mu)>"+str(ptCut)+" GeV, |#eta(#mu)|<"+str(etaCut)+"}")
text2.SetNDC()
text2.SetTextSize(0.04)
text2.Draw()
c3.SaveAs(args.saveDir+"summaryZStability"+suffix+".png")
c3.Close()

### corrected Z Counts per fill

graph_zcount=ROOT.TGraph(len(metaFills),metaFills,zcountl)
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

text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
text.SetNDC()
text.Draw()
text2=ROOT.TLatex(0.6,0.23,"#splitline{66 GeV<M(#mu#mu) < 116 GeV}{p_{T}(#mu)>"+str(ptCut)+" GeV, |#eta(#mu)|<"+str(etaCut)+"}")
text2.SetNDC()
text2.SetTextSize(0.04)
text2.Draw()
text3=ROOT.TLatex(0.6,0.33,"#color[4]{"+trigger+"}")
text3.SetNDC()
text3.SetTextSize(0.04)
text3.Draw()
c5.SaveAs(args.saveDir+"ZCountPerFill"+suffix+".png")
c5.Close()



### Accumalated corrected Z counts 

graph_zcountA=ROOT.TGraph(len(metaFills),timel,metazcountsAccu)
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
text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
text.SetNDC()
text.Draw()
text2=ROOT.TLatex(0.6,0.23,"#splitline{66 GeV<M(#mu#mu) < 116 GeV}{p_{T}(#mu)>"+str(ptCut)+" GeV, |#eta(#mu)|<"+str(etaCut)+"}")
text2.SetNDC()
text2.SetTextSize(0.04)
text2.Draw()
text3=ROOT.TLatex(0.6,0.33,"#color[4]{"+trigger+"}")
text3.SetNDC()
text3.SetTextSize(0.04)
text3.Draw()
c6.SaveAs(args.saveDir+"ZCountAccumulated"+suffix+".png")
c6.Close()
