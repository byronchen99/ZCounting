import ROOT
import pandas
import argparse

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cms", default="/eos/home-d/dwalter/www/ZCounting/CMS-2018-ZRateData/csvFiles/Mergedcsvfile.csv", type=str, help="give the CMS csv as input")
parser.add_argument("-s", "--saveDir", default='/eos/home-d/dwalter/www/ZCounting/CMS-2018-ZRateData/ZCrossSectionMonitoring/', type=str, help="give output dir")
args = parser.parse_args()

cmsfile=args.cms
# key words for columns of csv file
ZRate='ZRate'
Lumi='instDelLumi'

data = pandas.read_csv(cmsfile, sep=',')[[ZRate,Lumi]]
data['sigma'] = data[ZRate]/data[Lumi]

#remove outliers
data = data[data['sigma'] > 500]
data = data[data['sigma'] < 900]


graphXsecL=ROOT.TGraph(len(data),data[Lumi].values,data['sigma'].values)
graphXsecL.SetName("graph_metaXsecAtlas")
graphXsecL.SetMarkerStyle(23)
graphXsecL.SetMarkerColor(ROOT.kAzure-4)
graphXsecL.SetMarkerSize(1.5)
graphXsecL.SetTitle("Z cross section VS Lumi")
graphXsecL.GetXaxis().SetTitle("instantaneous luminosity [Hz/pb]")
graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")

print "the simple average cross section is "+str(sum(data['sigma'].values)/len(data))
graphXsecL.Fit("pol1","","",0.0045,0.018)
c3=ROOT.TCanvas("c3","c3",1000,600)
c3.SetGrid()
graphXsecL.Draw("AP")
c3.SaveAs(args.saveDir+""+ZRate+"_vs_"+Lumi+".png")
