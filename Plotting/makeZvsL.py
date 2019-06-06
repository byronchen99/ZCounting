import ROOT
import pandas
import numpy as np
import argparse
import pdb

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cms", default="/eos/home-d/dwalter/www/ZCounting/CMS-2018-ZRateData/csvFiles/Mergedcsvfile.csv", type=str, help="give the CMS csv as input")
parser.add_argument("-s", "--saveDir", default='/eos/home-d/dwalter/www/ZCounting/CMS-2018-ZRateData/ZCrossSectionMonitoring/', type=str, help="give output dir")
args = parser.parse_args()

cmsfile=args.cms
# key words for columns of csv file
#ZRate='ZRateUncorrected'
ZRate='ZRate'
Lumi='instDelLumi'

data = pandas.read_csv(cmsfile, sep=',')[['fill',ZRate,Lumi,'delZCount']]
data['sigma'] = data[ZRate]/data[Lumi]
data['sigmaE'] = data['sigma']/np.sqrt(data['delZCount']) 


# For each Run
for fill in data.drop_duplicates('fill')['fill'].values:
    dFill = data.loc[data['fill'] == fill]

    #pdb.set_trace()

    graphXsecL=ROOT.TGraphErrors(len(dFill),dFill[Lumi].values,dFill['sigma'].values,np.zeros(len(dFill)),dFill['sigmaE'].values)
    graphXsecL.SetName("graph_metaXsecAtlas")
    graphXsecL.SetMarkerStyle(23)
    graphXsecL.SetMarkerColor(ROOT.kAzure-4)
    graphXsecL.SetMarkerSize(1.5)
    graphXsecL.SetTitle("Z cross section VS Lumi")
    graphXsecL.GetXaxis().SetTitle("instantaneous luminosity [Hz/pb]")
    graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
    graphXsecL.GetYaxis().SetTitleOffset(1.0)

    print(">>> cross sections for Fill "+str(fill))
    print(">>> the simple average cross section is "+str(sum(data['sigma'].values)/len(data)))
    graphXsecL.Fit("pol1","","",0.0,0.025)
    c3=ROOT.TCanvas("c3","c3",1000,600)
    c3.SetGrid()
    graphXsecL.Draw("AP")

    legend=ROOT.TLegend(0.2,0.8,0.4,0.9)
    legend.AddEntry(graphXsecL,"Measurement (#pm stat.)","pe")
    legend.Draw("same")

    c3.SaveAs(args.saveDir+"/PlotsFill_"+str(fill)+"/"+ZRate+"_vs_"+Lumi+".png")


#remove outliers
data = data[data['sigma'] > 500]
data = data[data['sigma'] < 900]


graphXsecL=ROOT.TGraphErrors(len(data),data[Lumi].values,data['sigma'].values,np.zeros(len(data)),data['sigmaE'].values)
graphXsecL.SetName("graph_metaXsecAtlas")
graphXsecL.SetMarkerStyle(23)
graphXsecL.SetMarkerColor(ROOT.kAzure-4)
graphXsecL.SetMarkerSize(1.5)
graphXsecL.SetTitle("Z cross section VS Lumi")
graphXsecL.GetXaxis().SetTitle("instantaneous luminosity [Hz/pb]")
graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
graphXsecL.GetYaxis().SetTitleOffset(1.0)

print(">>> Producing all cross sections")
print(">>> the simple average cross section is "+str(sum(data['sigma'].values)/len(data)))
graphXsecL.Fit("pol1","","",0.0,0.025)
c3=ROOT.TCanvas("c3","c3",1000,600)
c3.SetGrid()
graphXsecL.Draw("AP")

legend=ROOT.TLegend(0.2,0.8,0.4,0.9)
legend.AddEntry(graphXsecL,"Measurement (#pm stat.)","pe")
legend.Draw("same")

c3.SaveAs(args.saveDir+""+ZRate+"_vs_"+Lumi+".png")
