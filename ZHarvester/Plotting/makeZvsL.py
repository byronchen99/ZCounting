import ROOT
import pandas as pd
import numpy as np
import argparse
import pdb
import os

pd.options.mode.chained_assignment = None

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)


latex = ROOT.TLatex()
latex.SetNDC()

parser = argparse.ArgumentParser()

parser.add_argument("cms", default="/eos/home-d/dwalter/www/ZCounting/CMS-2018-ZRateData/csvFiles/Mergedcsvfile.csv", type=str, help="give the CMS csv as input")
parser.add_argument("-s", "--saveDir", default='/eos/home-d/dwalter/www/ZCounting/CMS-2018-ZRateData/ZCrossSectionMonitoring/', type=str, help="give output dir")
args = parser.parse_args()

cmsfile=args.cms

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

ref_xsec = 610.1401700042975
ref_xsec_stat = 1.959
ref_xsec_lumi = 10.44

textsize = 0.03

# ------------------------------------------------------------------------------
def make_plots(ZXSec, xAxis='lumi',run_range=None, xTitle="instantaneous recorded luminosity [Hz/pb]", title="", plot_fills=False, nround=0):
    """
    valid xAxis: 'lumi', 'Pileup', 'measurement', 'time'
    """

    print(">>> load csv file in dataframe")
    data = pd.read_csv(cmsfile, sep=',')[['fill',ZXSec,'lumiRec','timewindow','pileUp','z_relstat','run','tdate_begin','tdate_end']]

    if run_range:
        data = data.query("run > {0} & run < {1}".format(*run_range))

    if xAxis == 'lumi':
        xAxis = 'lumiRecInst'
        data[xAxis] = data['lumiRec'] / data['timewindow']
    if xAxis == 'time':
        data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2
        #data['time'] -= 7200
        data['time'] -= (ROOT.TDatime(2042,01,01,12,00,00).Convert() - ROOT.TDatime(2017,01,01,12,00,00).Convert())


    print(">>> sort out points with low statistic")
    data = data.query("z_relstat < 0.05")

    data['sigma_E'] = data[ZXSec] * data['z_relstat']

    if plot_fills:
        print(">>> make plots for each fill")
        # For each Fill
        for fill in data.drop_duplicates('fill')['fill'].values:
            dFill = data.loc[data['fill'] == fill]

            subDir = outDir+"/PlotsXSec_"+str(fill)
            if not os.path.isdir(subDir):
                os.mkdir(subDir)

            graphXsecL=ROOT.TGraphErrors(len(dFill),dFill[xAxis].values,dFill[ZXSec].values,np.zeros(len(dFill)),data['sigma_E'].values)
            graphXsecL.SetName("graph_metaXsecAtlas")
            graphXsecL.SetMarkerStyle(23)
            graphXsecL.SetMarkerColor(ROOT.kAzure-4)
            graphXsecL.SetTitle(title)
            graphXsecL.GetXaxis().SetTitle(xTitle)
            graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
            graphXsecL.GetYaxis().SetTitleOffset(1.0)

            print(">>> cross sections for Fill "+str(fill))
            print(">>> the simple average cross section is "+str(sum(dFill[ZXSec].values)/len(dFill)))
            graphXsecL.Fit("pol1","","",0.0,0.025)
            c3=ROOT.TCanvas("c3","c3",1000,600)
            c3.SetGrid()
            graphXsecL.Draw("AP")


            legend=ROOT.TLegend(0.2,0.8,0.4,0.9)
            legend.AddEntry(graphXsecL,"Measurement (#pm stat.)","pe")
            legend.Draw("same")

            c3.SaveAs(subDir+"/"+ZXSec+"_vs_"+xAxis+".png")
            c3.Delete()


    print(">>> make plot for full data")
    if xAxis == 'measurement':
        xPoints = np.arange(0,len(data),1.)
    else:
        xPoints = data[xAxis].values

    graphXsecL=ROOT.TGraphErrors(len(data), xPoints, data[ZXSec].values, np.zeros(len(data)), data['sigma_E'].values)
    graphXsecL.SetName("graph_metaXsecAtlas")
    graphXsecL.SetMarkerStyle(33)
    graphXsecL.SetMarkerColor(ROOT.kAzure-4)
    graphXsecL.SetMarkerSize(1.)
    graphXsecL.SetTitle(title)
    graphXsecL.GetXaxis().SetTitle(xTitle)
    graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
    graphXsecL.GetYaxis().SetTitleOffset(1.0)



    #graphXsecL.GetYaxis().SetRangeUser(500.,700.)

    print(">>> Producing all cross sections")
    print(">>> the simple average cross section is "+str(sum(data[ZXSec].values)/len(data)))
    c3=ROOT.TCanvas("c3","c3",1000,600)
    c3.SetGrid()

    xmin = graphXsecL.GetXaxis().GetXmin()
    xmax = graphXsecL.GetXaxis().GetXmax()

    ref_line = ROOT.TLine(xmin,ref_xsec,xmax,ref_xsec)
    ref_line.SetLineColor(ROOT.kRed)
    ref_line.SetLineWidth(2)

    ref_box_inner = ROOT.TBox(xmin,ref_xsec-ref_xsec_stat,xmax,ref_xsec+ref_xsec_stat)
    ref_box_inner.SetFillColorAlpha(ROOT.kRed, 0.4)

    ref_box_outer = ROOT.TBox(xmin,ref_xsec-ref_xsec_lumi,xmax,ref_xsec+ref_xsec_lumi)
    ref_box_outer.SetFillColorAlpha(ROOT.kRed, 0.2)

    fit = graphXsecL.Fit("pol1","0S","",xmin,xmax)
    fit_a = fit.Value(0)
    fit_a_err = fit.ParError(0)
    fit_b = fit.Value(1)
    fit_b_err = fit.ParError(1)

    # line with statistical error
    xx = np.linspace(xmin,xmax,50)
    yy = fit_a + xx * fit_b
    yy_err_up = fit_a+fit_a_err + xx * (fit_b+fit_b_err)
    yy_err_down = fit_a-fit_a_err + xx * (fit_b-fit_b_err)
    fit_line = ROOT.TGraphAsymmErrors(len(xx), xx, yy, np.zeros(len(xx)), np.zeros(len(xx)), yy_err_down/yy, yy_err_up/yy)
    fit_line.SetFillColorAlpha(ROOT.kAzure-4,0.4)

    # line with error from luminosity linearity
    #lum_lin_err = 0.003 * 2 * abs(fit_b)
    #yy_err = xx * lum_lin_err
    #fit_line_lin = ROOT.TGraphErrors(len(xx), xx, yy, np.zeros(len(xx)), yy_err)
    #fit_line_lin.SetFillColor(ROOT.kAzure-4)

    if xAxis == 'time':
        graphXsecL.GetXaxis().SetTimeDisplay(1)
        fit_line.GetXaxis().SetTimeDisplay(1)

    entry = ROOT.TLegendEntry()
    entry.SetLabel("2017 High PU (#pm stat.))")
    entry.SetOption("pel")
    entry.SetTextSize(textsize)
    entry.SetTextFont(42)
    entry.SetMarkerStyle(33)
    entry.SetMarkerSize(1)
    entry.SetMarkerColor(ROOT.kAzure-4)
    entry.SetLineColor(ROOT.kBlack)
    entry.SetLineWidth(2)

    entry.SetTextAlign(12)

    entry2 = ROOT.TLegendEntry()
    entry2.SetLabel("2017 Low PU (#pm stat.) (#pm lumi.)")
    entry2.SetOption("fl")
    entry2.SetTextSize(0.03)
    entry2.SetTextFont(42)
    entry2.SetLineColor(ROOT.kRed)
    #entry2.SetLineWidth(2)
    entry2.SetFillColorAlpha(ROOT.kRed, 0.4)
    entry2.SetTextAlign(12)
    entry2.SetMarkerSize(0)

    legend=ROOT.TLegend(0.2,0.8,0.5,0.9)
    legend.AddEntry(graphXsecL, "2017 High PU (#pm stat.)", "pe")
    #legend.GetListOfPrimitives().Add(entry)
    legend.GetListOfPrimitives().Add(entry2)
    legend.SetTextSize(textsize)
    legend.SetTextFont(42)

    # --- Plot

    graphXsecL.Draw("AP")

    ref_line.Draw("same")
    ref_box_inner.Draw("same")
    ref_box_outer.Draw("same")

    fit_line.Draw("3 same")
    #fit_line_lin.Draw("3 same")

    graphXsecL.Draw("PEsame")


    latex.SetTextAlign(11)
    latex.SetTextFont(42)
    latex.SetTextSize(textsize)
    latex.DrawLatex(0.6, 0.89, "Fit: y = a + b \cdot x")
    latex.DrawLatex(0.6, 0.85, "a =\ {0} \pm {1} (stat.)".format(round(fit_a,1), round(fit_a_err,1)))
    if nround:
        latex.DrawLatex(0.6, 0.81, "b = {0} \pm {1} (stat.)".format(round(fit_b,nround), round(fit_b_err,nround)))
        #latex.DrawLatex(0.6, 0.81, "b = {0} \pm {1} (stat.) \pm {2} (lumi.)".format(round(fit_b,nround), round(fit_b_err,nround), round(lum_lin_err,nround)))
    else:
        latex.DrawLatex(0.6, 0.81, "b = {0} \pm {1} (stat.)".format(int(round(fit_b)), int(round(fit_b_err))))
        #latex.DrawLatex(0.6, 0.81, "b = {0} \pm {1} (stat.) \pm {2} (lumi.)".format(int(round(fit_b)), int(round(fit_b_err)), int(round(lum_lin_err))))

    latex.DrawLatex(0.6, 0.77, "\chi^2 / ndf = {0}".format(round(fit.Chi2() / len(data),2)))


    legend.Draw("same")

    outstring = "{0}/{1}_vs_{2}".format(outDir,ZXSec,xAxis)
    if run_range:
        outstring += "_run{0}to{1}".format(*run_range)

    c3.SaveAs(outstring+".png")
    c3.Close()

# ------------------------------------------------------------------------------
make_plots("zXSec", title="Z rate / lumi - 2017B-F")#, plot_fills=True)

make_plots("zXSec_mc", title="Corrected Z rate / lumi - 2017B-F", plot_fills=True)
make_plots("zXSec_mc", title="Corrected Z rate / lumi - 2017B-E", run_range=(297046,304797))
make_plots("zXSec_mc", title="Corrected Z rate / lumi - 2017B", run_range=(297046,299329))
make_plots("zXSec_mc", title="Corrected Z rate / lumi - 2017C", run_range=(299368,302029))
make_plots("zXSec_mc", title="Corrected Z rate / lumi - 2017D", run_range=(302030,303434))
make_plots("zXSec_mc", title="Corrected Z rate / lumi - 2017E", run_range=(303434,304797))
make_plots("zXSec_mc", title="Corrected Z rate / lumi - 2017F", run_range=(305040,306462))

make_plots("zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi", nround=4)
make_plots("zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017B", run_range=(297046,299329), nround=4)
make_plots("zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017C", run_range=(299368,302029), nround=4)
make_plots("zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017D", run_range=(302030,303434), nround=4)
make_plots("zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017E", run_range=(303434,304797), nround=4)
make_plots("zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017F", run_range=(305040,306462), nround=4)
