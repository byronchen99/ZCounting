import os,sys
#sys.argv.append( '-b-' )
#from scipy.interpolate import UnivariateSpline
import ROOT
from ROOT import TGraph
from ROOT import TGraphErrors
from array import array
from ROOT import *
from operator import truediv
import random
import argparse
import scipy.integrate as integrate
import pandas
import csv
import math
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
#ROOT.gStyle.SetTitleAlign(33)
ROOT.gStyle.SetTitleX(.3)

ZeffUnc=0.03

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--atlas", default="nothing", type=string, help="give the ATLAS csv as input")
parser.add_argument("-c", "--cms", default="nothing", type=string, help="give the CMS csv as input")
parser.add_argument("-f", "--fill", default='', help="give fill numbers")
parser.add_argument("-s", "--saveDir", default='', help="give save dir")
args = parser.parse_args()

if args.cms=="nothing":
	print "please provide cms input file"
	sys.exit()
if args.atlas=="nothing":
	print "please provide atlas input file"
	sys.exit()

print args.cms
print args.atlas

atlasfile=open(str(args.atlas))
cmsfile=open(str(args.cms))
saveDir=args.saveDir
fills=args.fill.split(",")
afills=[]

if fills==['']:
	data=pandas.read_csv(str(args.cms), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
	print data.axes
	fills=data.drop_duplicates('fill')['fill'].tolist()

	data2=pandas.read_csv(str(args.atlas), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
	print data2.axes
	afills=data2.drop_duplicates('fill')['fill'].tolist()

print "Fills being processed: "
print fills

suffix=""
if "Central" in str(args.cms):
	suffix="Barrel"
else:
	suffix="Inclusive"

linescms=cmsfile.readlines()
linesatlas=atlasfile.readlines()
metaFills=array('d')
metaXsecCMS=array('d')
metaXsecATLAS=array('d')
metaZLumiRatio=array('d')
metaZLumiRatioEx=array('d')
metaZLumiRatioEy=array('d')
metaLines = [["Fill","Z-Ratio","Z-Ratio uncertainty"]]

for fill in fills:
	if fill not in afills:
		print str(fill)+"  not in atlas file"
		continue
#	print fill
	cmsRates=array('d')
	cmsRatesE=array('d')
	cmsTimes=array('d')
	cmsTimesE=array('d')
	cmsInstLum=array('d')
	cmsXsec=array('d')
	cmsXsecEx=array('d')
	cmsXsecEy=array('d')

	atlasRates=array('d')
	atlasRatesE=array('d')
	atlasTimes=array('d')
	atlasTimesE=array('d')
	atlasInstLum=array('d')
	atlasXsec=array('d')
	atlasXsecEx=array('d')
	atlasXsecEy=array('d')

	ZintRatio=array('d')

	atlascmsratio=array('d')
	atlascmsratioerrorX=array('d')
	atlascmsratioerrorXE=array('d')
	atlascmsratioerrorY=array('d')
	atlascmsratioerrorYE=array('d')
	
	k=0
	for linecms in range(0,len(linescms)):
		elements=linescms[linecms].split(",")
		
		if elements[0]==str(fill):
			if float(elements[6])<1000.:
				continue
			k=k+1
			rate=elements[3]
			cmsRates.append(float(rate)) #*1.13
			cmsRatesE.append( float(rate) * math.sqrt( ( math.sqrt(float(elements[6]))/float(elements[6]) )**2 + ZeffUnc**2 ) )
			datestamp=elements[1].split(" ")
			date=ROOT.TDatime(2018,int(datestamp[0].split("/")[0]),int(datestamp[0].split("/")[1]),int(datestamp[1].split(":")[0]),int(datestamp[1].split(":")[1]),int(datestamp[1].split(":")[2]))
			cmsTimes.append(date.Convert())
			cmsTimesE.append(date.Convert()-date.Convert())
			cmsInstLum.append(float(elements[4]))
			cmsXsec.append(float(elements[3])/float(elements[4]))
			#cmsXsecEy.append((float(elements[6])/float(elements[5]))*0.02)
			cmsXsecEy.append(math.sqrt( ( math.sqrt(float(elements[6]))/float(elements[6]) )**2 + ZeffUnc**2 ))
			
	if len(cmsTimes)==0:
		continue

	graph_cms=ROOT.TGraphErrors(k,cmsTimes,cmsRates,cmsTimesE,cmsRatesE)
	graph_cms.SetName("graph_cms")
	graph_cms.SetMarkerStyle(22)
	graph_cms.SetMarkerColor(kOrange+8)
	graph_cms.SetFillStyle(0)
	graph_cms.SetMarkerSize(1.5)

	graph_cmsinstLum=ROOT.TGraph(k,cmsTimes,cmsInstLum)
	graph_cmsinstLum.SetName("graph_cmsinstLum")
	graph_cmsinstLum.SetMarkerStyle(34)
	graph_cmsinstLum.SetMarkerColor(kRed)
	graph_cmsinstLum.SetMarkerSize(1.5)

	graph_cmsXsec=ROOT.TGraph(k,cmsTimes,cmsXsec)
	graph_cmsXsec.SetName("graph_cmsXsec")
	graph_cmsXsec.SetMarkerStyle(22)
	graph_cmsXsec.SetMarkerColor(kOrange+8)
	graph_cmsXsec.SetFillStyle(0)
	graph_cmsXsec.SetMarkerSize(1.5)

	i=0
	for lineatlas in range(0,len(linesatlas)-1):
		elements=linesatlas[lineatlas].split(",")
		if elements[0]==str(fill):
			i=i+1
			rate=elements[3]
			atlasRates.append(float(rate))
			atlasRatesE.append( float(rate) * math.sqrt( ( math.sqrt(float(elements[6]))/float(elements[6]) )**2 + ZeffUnc**2 )  )
			datestamp=elements[1].split(" ")
			date=ROOT.TDatime(2018,int(datestamp[1].split("/")[1]),int(datestamp[1].split("/")[2]),int(datestamp[2].split(":")[0]),int(datestamp[2].split(":")[1]),int(datestamp[2].split(":")[2]))
			atlasTimes.append(date.Convert())
			atlasTimesE.append(date.Convert()-date.Convert())
			atlasInstLum.append(float(elements[4]))
			if graph_cms.Eval(date.Convert())!=0:
				atlascmsratio.append(float(rate)/graph_cms.Eval(date.Convert()))
				atlascmsratioerrorY.append(float(rate)/graph_cms.Eval(date.Convert()))
			else:
				atlascmsratio.append(1.2) #some value	
				atlascmsratioerrorY.append(1.2) #some value	
			atlasXsec.append(float(elements[3])/float(elements[4]))
			atlasXsecEy.append(math.sqrt( ( math.sqrt(float(elements[6]))/float(elements[6]) )**2 + ZeffUnc**2 ))
			atlascmsratioerrorYE.append(0.04)
			atlascmsratioerrorX.append(date.Convert())	
			atlascmsratioerrorXE.append(0.0)	

	if len(atlasTimes)==0:
		continue		
	
	graph_atlas=ROOT.TGraphErrors(i,atlasTimes,atlasRates,atlasTimesE,atlasRatesE)
	graph_atlas.SetName("graph_atlas")
	graph_atlas.SetMarkerStyle(23)
	graph_atlas.SetFillStyle(0)
	graph_atlas.SetMarkerColor(kAzure-4)
	graph_atlas.SetMarkerSize(1.5)
	
	graph_ratioE=ROOT.TGraphErrors(i,atlascmsratioerrorX,atlascmsratioerrorY,atlascmsratioerrorXE,atlascmsratioerrorYE)
	graph_ratioE.SetName("graph_ratioE")
	graph_ratioE.SetFillColor(kOrange)
	graph_ratioE.SetFillStyle(3001)
	graph_ratioE.SetTitle("")
	graph_ratioE.SetMarkerStyle(20)
	graph_ratioE.SetMarkerColor(kGray+1)
	graph_ratioE.SetMarkerSize(1.5)
	graph_ratioE.GetXaxis().SetTimeDisplay(1)
	graph_ratioE.GetXaxis().SetLabelSize(0.1)
	graph_ratioE.GetYaxis().SetRangeUser(0.75,1.25)
	graph_ratioE.GetYaxis().SetTitle("ATLAS/CMS")
	graph_ratioE.GetYaxis().SetTitleSize(0.1)
	graph_ratioE.GetYaxis().SetTitleOffset(0.45)
	graph_ratioE.GetYaxis().SetLabelSize(0.08)

	graph_atlascmsratio=ROOT.TGraph(i,atlasTimes,atlascmsratio)
	graph_atlascmsratio.SetName("graph_atlascmsratio")
	graph_atlascmsratio.SetTitle("")
	graph_atlascmsratio.SetMarkerStyle(20)
	graph_atlascmsratio.SetMarkerColor(kGray+1)
	graph_atlascmsratio.SetMarkerSize(1.5)
	graph_atlascmsratio.GetXaxis().SetTimeDisplay(1)
	graph_atlascmsratio.GetXaxis().SetLabelSize(0.1)
	graph_atlascmsratio.GetYaxis().SetRangeUser(0.75,1.25)
	graph_atlascmsratio.GetYaxis().SetTitle("ATLAS/CMS")
	graph_atlascmsratio.GetYaxis().SetTitleSize(0.1)
	graph_atlascmsratio.GetYaxis().SetTitleOffset(0.45)
	graph_atlascmsratio.GetYaxis().SetLabelSize(0.08)

	graph_atlasinstLum=ROOT.TGraph(i,atlasTimes,atlasInstLum)
	graph_atlasinstLum.SetName("graph_atlasinstLum")
	graph_atlasinstLum.SetMarkerStyle(34)
	graph_atlasinstLum.SetMarkerColor(kBlue)
	graph_atlasinstLum.SetMarkerSize(2)

	def atlasFunc(x):
		return graph_atlas.Eval(x)

	def cmsFunc(x):
		return graph_cms.Eval(x)


	startTime=cmsTimes[0]
	endTime=cmsTimes[-1]
	boundUpTime=atlasTimes[-1]
	boundDownTime=atlasTimes[0]
	if atlasTimes[-1]<cmsTimes[-1]:
		endTime=atlasTimes[-1]
		boundUpTime=cmsTimes[-1]
	if atlasTimes[0]>cmsTimes[0]:
		startTime=atlasTimes[0]
		boundDownTime=cmsTimes[0]

	cmsTimesZInt=array('d')
	for cmsTime in range(0,len(cmsTimes)):
		if cmsTimes[cmsTime]>startTime and cmsTimes[cmsTime]<=endTime:
			ZintRatio.append(integrate.quad(atlasFunc,startTime,cmsTimes[cmsTime])[0]/integrate.quad(cmsFunc,startTime,cmsTimes[cmsTime])[0])	
			cmsTimesZInt.append(cmsTimes[cmsTime])
			#ZintRatio.append(integrate.quad(atlasFunc,cmsTimes[0],cmsTimes[cmsTime])[0]/integrate.quad(cmsFunc,cmsTimes[0],cmsTimes[cmsTime])[0])					
	#print "CMS Error: "+str(integrate.quad(atlasFunc,cmsTimes[0],cmsTimes[cmsTime])[1]/integrate.quad(atlasFunc,cmsTimes[0],cmsTimes[cmsTime])[0])	
	#print "ATLAS Error: "+str(integrate.quad(cmsFunc,cmsTimes[0],cmsTimes[cmsTime])[1]/integrate.quad(cmsFunc,cmsTimes[0],cmsTimes[cmsTime])[0])	

	if len(ZintRatio)<1:
		continue
	
	graph_ZintRatio=ROOT.TGraph(len(cmsTimesZInt),cmsTimesZInt,ZintRatio)
	graph_ZintRatio.SetName("graph_ZintRatio")
	graph_ZintRatio.SetMarkerStyle(34)
	
	graph_ZintRatio.SetLineColor(kBlue)
	graph_ZintRatio.SetLineWidth(3)
	graph_ZintRatio.SetMarkerSize(2)

	c1=ROOT.TCanvas("c1","c1",1000,600)
	c1.Divide(1,2)
	c1.cd(1).SetPad(0.0,0.3,1.0,1.0)

	graph_atlas.GetXaxis().SetTimeDisplay(1)
	graph_atlas.SetTitle(suffix+" Z-Rates, Fill "+str(fill))
	graph_atlas.GetYaxis().SetTitle("Z-Rate [Hz]")
	graph_atlas.GetYaxis().SetTitleSize(0.07)
	graph_atlas.GetYaxis().SetTitleOffset(0.5)
	graph_atlas.GetXaxis().SetTitle("Time")
	graph_atlas.GetXaxis().SetTitleSize(0.06)
	graph_atlas.GetXaxis().SetTitleOffset(0.72)
	graph_atlas.GetXaxis().SetLabelSize(0.05)
	graph_atlas.GetYaxis().SetLabelSize(0.05)
	graph_atlas.GetXaxis().SetRangeUser(boundDownTime,boundUpTime)
	
	graph_cms.GetXaxis().SetTimeDisplay(1)
	graph_cms.SetTitle(suffix+" Z-Rates, Fill "+str(fill))
	graph_cms.GetYaxis().SetTitle("Z-Rate [Hz]")
	graph_cms.GetYaxis().SetTitleSize(0.07)
	graph_cms.GetYaxis().SetTitleOffset(0.5)
	graph_cms.GetXaxis().SetTitle("Time")
	graph_cms.GetXaxis().SetTitleSize(0.06)
	graph_cms.GetXaxis().SetTitleOffset(0.72)
	graph_cms.GetXaxis().SetLabelSize(0.05)
	graph_cms.GetYaxis().SetLabelSize(0.05)
	graph_cms.GetXaxis().SetRangeUser(boundDownTime,boundUpTime)
		
	c1.cd(1)
	graph_atlas.Draw("AP")
	graph_cms.Draw("Psame")	

	legend=ROOT.TLegend(0.65,0.65,0.9,0.9)
	legend.AddEntry(graph_cms,"CMS","p")
	legend.AddEntry(graph_atlas,"ATLAS","p")
	legend.AddEntry(graph_ratioE,"Ratio","AP3Elf")
	legend.AddEntry(graph_ZintRatio,"Int. Z-Luminosity Ratio","l")
	legend.Draw()

	text1=ROOT.TText(0.2,0.85,"Work In Progress")
	text1.SetNDC()
	text1.Draw()

	c1.cd(2).SetPad(0.0,0.0,1.0,0.3)
	c1.cd(2)
	text2=ROOT.TText(0.8,0.73,"Int. Z-Luminosity Ratio "+str(ZintRatio[-1]))
	text2.SetNDC()
	text2.Draw()
		
	graph_ratioE.Draw("AP3")	
	graph_ratioE.GetXaxis().SetRangeUser(boundDownTime,boundUpTime)
	graph_atlascmsratio.GetXaxis().SetRangeUser(boundDownTime,boundUpTime)

	line=ROOT.TLine(boundDownTime,1,boundUpTime,1)
	lineUp=ROOT.TLine(boundDownTime,1.1,boundUpTime,1.1)
	lineDown=ROOT.TLine(boundDownTime,0.9,boundUpTime,0.9)

	lineUp.SetLineStyle(2)
	lineDown.SetLineStyle(2)
	line.Draw("same")
	lineUp.Draw("same")
	lineDown.Draw("same")
	graph_ZintRatio.Draw("samel")
	graph_ZintRatio.GetXaxis().SetRangeUser(boundDownTime,boundUpTime)

        text2=ROOT.TText(0.67,0.79,"Int. Z-Luminosity Ratio "+str(round(ZintRatio[-1],3)))
	text2.SetNDC()
	text2.SetTextSize(0.1)
	text2.Draw()

	c1.cd(1)
	c1.Update()
	#c1.SaveAs("zrates"+str(fill)+suffix+".root")
	c1.SaveAs(saveDir+"/ATLAS_CMS_Comparison/zrates"+str(fill)+suffix+".png")
	c1.Close()
	
	metaXsecCMS.append(sum(cmsXsec)/len(cmsXsec))
	metaXsecATLAS.append(sum(atlasXsec)/len(atlasXsec))
	metaZLumiRatio.append(ZintRatio[-1])
	metaZLumiRatioEx.append(0)
	metaZLumiRatioEy.append(ZintRatio[-1]*math.sqrt(ZeffUnc**2 + ZeffUnc**2))
	metaFills.append(float(fill))	

	metaLines.append([fill,ZintRatio[-1],ZintRatio[-1]*math.sqrt(ZeffUnc**2 + ZeffUnc**2)])

	atlasXsec2=array('d')
	for n in range(0,len(atlasXsec)):
		atlasXsec2.append(atlasXsec[n]/(sum(atlasXsec)/len(atlasXsec)))

	graph_atlasXsec2=ROOT.TGraphErrors(i,atlasTimes,atlasXsec2,atlasTimesE,atlasXsecEy)
	graph_atlasXsec2.SetName("graph_atlasXsec2")
	graph_atlasXsec2.SetTitle(suffix+" Z-Rates, Fill "+str(fill))
	graph_atlasXsec2.SetMarkerStyle(23)
	graph_atlasXsec2.SetFillStyle(0)
	graph_atlasXsec2.SetMarkerColor(kAzure-4)
	graph_atlasXsec2.SetMarkerSize(2)
	graph_atlasXsec2.GetXaxis().SetTimeDisplay(1)
	graph_atlasXsec2.GetYaxis().SetTitle("#sigma^{fid}_{Z}/<#sigma^{fid}_{Z}>")
	graph_atlasXsec2.GetYaxis().SetTitleSize(0.05)
	graph_atlasXsec2.GetYaxis().SetTitleOffset(0.95)
	graph_atlasXsec2.GetXaxis().SetTitle("Time")
	graph_atlasXsec2.GetXaxis().SetTitleSize(0.06)
	graph_atlasXsec2.GetXaxis().SetTitleOffset(0.72)
	graph_atlasXsec2.GetXaxis().SetLabelSize(0.05)
	graph_atlasXsec2.GetYaxis().SetLabelSize(0.05)	
	graph_atlasXsec2.GetYaxis().SetRangeUser(0.9,1.1)

	cmsXsec2=array('d')
	for n in range(0,len(cmsXsec)):
		cmsXsec2.append(cmsXsec[n]/(sum(cmsXsec)/len(cmsXsec)))	

	graph_cmsXsec2=ROOT.TGraphErrors(len(cmsXsec),cmsTimes,cmsXsec2,cmsTimesE,cmsXsecEy)
	graph_cmsXsec2.SetName("graph_cmsXsec2")
	graph_cmsXsec2.SetTitle(suffix+" Z-Rates, Fill "+str(fill))
	graph_cmsXsec2.SetMarkerStyle(22)
	graph_cmsXsec2.SetMarkerColor(kOrange+8)
	graph_cmsXsec2.SetFillStyle(0)
	graph_cmsXsec2.SetMarkerSize(2)
	graph_cmsXsec2.GetXaxis().SetTimeDisplay(1)
	graph_cmsXsec2.GetYaxis().SetTitle("#sigma^{fid}_{Z}/<#sigma^{fid}_{Z}>")
	graph_cmsXsec2.GetYaxis().SetTitleSize(0.05)
	graph_cmsXsec2.GetYaxis().SetTitleOffset(0.95)
	graph_cmsXsec2.GetXaxis().SetTitle("Time")
	graph_cmsXsec2.GetXaxis().SetTitleSize(0.06)
	graph_cmsXsec2.GetXaxis().SetTitleOffset(0.72)
	graph_cmsXsec2.GetXaxis().SetLabelSize(0.05)
	graph_cmsXsec2.GetYaxis().SetLabelSize(0.05)	
	graph_cmsXsec2.GetYaxis().SetRangeUser(0.9,1.1)
	
	c4=ROOT.TCanvas("c4","c4",1000,600)
	c4.SetGrid()
	graph_atlasXsec2.Draw("AP")
	graph_cmsXsec2.Draw("Psame")	
		
	legend=ROOT.TLegend(0.75,0.75,0.9,0.9)
	legend.AddEntry(graph_cmsXsec2,"CMS","p")
	legend.AddEntry(graph_atlasXsec2,"ATLAS","p")
	legend.Draw()
	text=ROOT.TText(0.2,0.85,"Work In Progress")
	text.SetNDC()
	text.Draw()
	#c4.SaveAs("ZStability"+str(fill)+suffix+".root")
	c4.SaveAs(saveDir+"/ATLAS_CMS_Comparison/ZStability"+str(fill)+suffix+".png")	
	c4.Close()
		
ROOT.gROOT.SetBatch(True)

metaXsecATLAS2=array('d')
for n in range(0,len(metaXsecATLAS)):
	metaXsecATLAS2.append(metaXsecATLAS[n]/(sum(metaXsecATLAS)/len(metaXsecATLAS)))	

metaXsecCMS2=array('d')
for n in range(0,len(metaXsecCMS)):
	metaXsecCMS2.append(metaXsecCMS[n]/(sum(metaXsecCMS)/len(metaXsecCMS)))	

graph_metaatlasXsec=ROOT.TGraph(len(metaFills),metaFills,metaXsecATLAS2)
graph_metaatlasXsec.SetName("graph_metaXsecAtlas")
graph_metaatlasXsec.SetMarkerStyle(23)
graph_metaatlasXsec.SetMarkerColor(kAzure-4)
graph_metaatlasXsec.SetMarkerSize(1.5)
graph_metaatlasXsec.SetTitle(suffix+" Z-Rates")

graph_metacmsXsec=ROOT.TGraph(len(metaFills),metaFills,metaXsecCMS2)
graph_metacmsXsec.SetName("graph_metaXsecCms")
graph_metacmsXsec.SetMarkerStyle(22)
graph_metacmsXsec.SetMarkerColor(kOrange+8)
graph_metacmsXsec.SetMarkerSize(1.5)
graph_metacmsXsec.SetTitle(suffix+" Z-Rates")

multMetaGraphXsec=ROOT.TMultiGraph("multMetaGraphXsec",suffix+" Z-Rates")
multMetaGraphXsec.SetName("multMetaGraphXsec")

graph_metacmsXsec.GetXaxis().SetTitle("Fill")
graph_metacmsXsec.GetYaxis().SetTitle("#sigma^{fid}_{Z}/<#sigma^{fid}_{Z}>")
graph_metacmsXsec.GetXaxis().SetTitleSize(0.05)
graph_metacmsXsec.GetYaxis().SetTitleSize(0.05)
graph_metacmsXsec.GetXaxis().SetTitleOffset(0.87)
graph_metacmsXsec.GetYaxis().SetTitleOffset(0.87)

graph_metaatlasXsec.GetXaxis().SetTitle("Fill")
graph_metaatlasXsec.GetYaxis().SetTitle("#sigma^{fid}_{Z}/<#sigma^{fid}_{Z}>")
graph_metaatlasXsec.GetXaxis().SetTitleSize(0.05)
graph_metaatlasXsec.GetYaxis().SetTitleSize(0.05)
graph_metaatlasXsec.GetXaxis().SetTitleOffset(0.87)
graph_metaatlasXsec.GetYaxis().SetTitleOffset(0.87)

multMetaGraphXsec.Add(graph_metacmsXsec)
multMetaGraphXsec.Add(graph_metaatlasXsec)

c3=ROOT.TCanvas("c3","c3",1000,600)
c3.SetGrid()

graph_metaatlasXsec.Draw("AP")
graph_metacmsXsec.Draw("Psame")

if suffix=="Barrel":
	graph_metaatlasXsec.GetYaxis().SetRangeUser(0.9,1.1)
if suffix=="Inclusive":
	graph_metaatlasXsec.GetYaxis().SetRangeUser(0.9,1.1)

legend=ROOT.TLegend(0.75,0.75,0.9,0.9)
legend.AddEntry(graph_metacmsXsec,"CMS","p")
legend.AddEntry(graph_metaatlasXsec,"ATLAS","p")
legend.Draw()

text=ROOT.TText(0.2,0.85,"Work In Progress")
text.SetNDC()
text.Draw()

#c3.SaveAs("summaryZStability"+suffix+".root")
c3.SaveAs(saveDir+"/summaryZStability"+suffix+"AtlasCMS.png")
c3.Close()

graph_metaZIntLumi=ROOT.TGraphErrors(len(metaFills),metaFills,metaZLumiRatio,metaZLumiRatioEx,metaZLumiRatioEy)
graph_metaZIntLumi.SetName("graph_metaZIntLumi")
graph_metaZIntLumi.SetTitle(suffix+" Z-Rates")
graph_metaZIntLumi.SetMarkerStyle(34)
graph_metaZIntLumi.SetMarkerColor(kBlue)
graph_metaZIntLumi.SetMarkerSize(2)
graph_metaZIntLumi.GetXaxis().SetTitle("Fill")
graph_metaZIntLumi.GetYaxis().SetTitle("Z Int. Luminosity Ratio ATLAS/CMS")
graph_metaZIntLumi.GetXaxis().SetTitleSize(0.05)
graph_metaZIntLumi.GetYaxis().SetTitleSize(0.05)
graph_metaZIntLumi.GetXaxis().SetTitleOffset(0.87)
graph_metaZIntLumi.GetYaxis().SetTitleOffset(0.87)

c5=ROOT.TCanvas("c5","c5",1000,600)
c5.SetGrid()

graph_metaZIntLumi.Draw("APE")
graph_metaZIntLumi.GetYaxis().SetRangeUser(0.85,1.15)
text=ROOT.TText(0.2,0.85,"Work In Progress")
text.SetNDC()
text.Draw()
#c5.SaveAs(suffix+"summaryZIntLumiRatio"+suffix+".root")
c5.SaveAs(saveDir+"/"+suffix+"summaryZIntLumiRatio"+suffix+".png")
c5.Close()

myFile = open(saveDir+'/ATLAS_CMS_ZRatiosPerFill_2018_v1.csv', 'w')
with myFile:
    writer = csv.writer(myFile)
    writer.writerows(metaLines)
     
print("Writing complete")
