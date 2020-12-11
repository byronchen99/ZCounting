import ROOT
import pandas as pd
import numpy as np
import argparse
import pdb
import os
import math

pd.options.mode.chained_assignment = None

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)


latex = ROOT.TLatex()
latex.SetNDC()

parser = argparse.ArgumentParser()

parser.add_argument("--rates", required=True, type=str, help="csv file with z rates per Measurement")
parser.add_argument("--xsec",  type=str, help="csv file where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")

args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

plot_lowPU = False
if args.xsec:
    plot_lowPU = True
    # --- get Z xsec
    data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
    xsecBB = sum(data_xsec['zDelBB_mc'])/sum(data_xsec['lumiRec'])
    xsecBE = sum(data_xsec['zDelBE_mc'])/sum(data_xsec['lumiRec'])
    xsecEE = sum(data_xsec['zDelEE_mc'])/sum(data_xsec['lumiRec'])


textsize = 0.03

# ------------------------------------------------------------------------------
def make_plots(data, ZXSec, relstat,
    xAxis='lumi',run_range=None,
    xTitle="instantaneous recorded luminosity [Hz/pb]",
    title="", year="2017", region='inclusive',
    plot_fills=False, plot_all=False,
    nround=0, normalized=False
):
    """
    valid xAxis: 'lumi', 'Pileup', 'measurement', 'time'
    """

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
    data = data.query(relstat+" < 0.05")

    if normalized:
        data[ZXSec] /= data[ZXSec].mean()

    data['sigma_E'] = data[ZXSec] * data[relstat]


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
            if normalized:
                graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} / <#sigma^{fid}_{Z}>")
            else:
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

    if plot_all:
        print(">>> make plot all measurements")
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
        if normalized:
            graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} / <#sigma^{fid}_{Z}>")
        else:
            graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
        graphXsecL.GetYaxis().SetTitleOffset(1.0)

        #graphXsecL.GetYaxis().SetRangeUser(500.,700.)

        print(">>> Producing all cross sections")
        print(">>> the simple average cross section is "+str(sum(data[ZXSec].values)/len(data)))
        c3=ROOT.TCanvas("c3","c3",1000,600)
        c3.SetGrid()

        xmin = graphXsecL.GetXaxis().GetXmin()
        xmax = graphXsecL.GetXaxis().GetXmax()

        if plot_lowPU:
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

        print("fit: ",fit_a, fit_b)

        # line with statistical error
        xx = np.linspace(xmin,xmax,50)
        yy = fit_a + xx * fit_b
        yy_err_up = fit_a+fit_a_err + xx * (fit_b+fit_b_err)
        yy_err_down = fit_a-fit_a_err + xx * (fit_b-fit_b_err)
        fit_line = ROOT.TGraphAsymmErrors(len(xx), xx, yy, np.zeros(len(xx)), np.zeros(len(xx)), yy - yy_err_down, yy_err_up - yy)
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
        entry.SetLabel("Measurement (#pm stat.))")
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
        if plot_lowPU:
            legend.GetListOfPrimitives().Add(entry2)
        legend.SetTextSize(textsize)
        legend.SetTextFont(42)

        # --- Plot

        graphXsecL.Draw("AP")

        if args.xsec:
            ref_line.Draw("same")
            ref_box_inner.Draw("same")
            ref_box_outer.Draw("same")

        fit_line.Draw("3 same")
        #fit_line_lin.Draw("3 same")

        graphXsecL.Draw("PEsame")


        latex.SetTextAlign(11)
        latex.SetTextFont(42)
        latex.SetTextSize(textsize)
        latex.DrawLatex(0.7, 0.89, "Fit: y = a + b \cdot x")

        from math import log10, floor
        def round_sig(x, sig=2):
            return round(x, sig-int(floor(log10(abs(x))))-1)


        latex.DrawLatex(0.7, 0.85, "a =\ {0} \pm {1} (stat.)".format(round_sig(fit_a), round_sig(fit_a_err)))
        latex.DrawLatex(0.7, 0.81, "b = {0} \pm {1} (stat.)".format(round_sig(fit_b), round_sig(fit_b_err)))

        latex.DrawLatex(0.7, 0.77, "\chi^2 / ndf = {0}".format(round_sig(fit.Chi2() / len(data),3)))


        legend.Draw("same")

        outstring = "{0}/{1}_vs_{2}".format(outDir,ZXSec,xAxis)
        if run_range:
            outstring += "_run{0}to{1}".format(*run_range)

        c3.SaveAs(outstring+".png")
        c3.Close()

    print(">>> make plot and combine measurements into bins")
    x_step = 0.0005
    xx = np.arange(data[xAxis].min()*(1./x_step)//1/(1./x_step),data[xAxis].max()*(1./x_step)//1/(1./x_step), x_step)
    xx = np.append(xx, xx[-1]+x_step)
    xx_centers = []
    yy = []
    yy_err = []
    for i in range(0,len(xx)-1):
        dyy = data.loc[(data[xAxis] < xx[i+1]) & (data[xAxis] > xx[i])]
        if len(dyy)==0:
            continue
        xx_centers.append((xx[i] + (xx[i+1] - xx[i] / 2)))
        yy_avg = np.average(dyy[ZXSec].values, weights=1./dyy['sigma_E'].values)
        yy.append(yy_avg)
        yy_err.append(np.sqrt(1./sum((1./dyy['sigma_E'].values)**2)))

    graphXsecL=ROOT.TGraphErrors(len(xx_centers), np.array(xx_centers), np.array(yy), np.ones(len(xx_centers))*x_step/2, np.array(yy_err))
    graphXsecL.SetName("graph_metaXsecAtlas")
    graphXsecL.SetMarkerStyle(33)
    graphXsecL.SetMarkerColor(ROOT.kAzure-4)
    graphXsecL.SetMarkerSize(1.)
    graphXsecL.SetTitle("")
    graphXsecL.GetXaxis().SetTitle(xTitle)
    if normalized:
        graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} / <#sigma^{fid}_{Z}>")
    else:
        graphXsecL.GetYaxis().SetTitle("#sigma^{fid}_{Z} [pb]")
    graphXsecL.GetYaxis().SetTitleOffset(1.2)

    xmin = min(xx_centers) - x_step
    xmax = max(xx_centers) + x_step
    ymin = min(np.array(yy)-np.array(yy_err))
    ymax = max(np.array(yy)+np.array(yy_err))
    graphXsecL.GetXaxis().SetRangeUser(xmin, xmax)
    graphXsecL.GetYaxis().SetRangeUser(ymin - 0.05*(ymax-ymin), ymax + 0.2*(ymax-ymin))

    fit = graphXsecL.Fit("pol1","0S","",xmin,xmax)
    fit_a = fit.Value(0)
    fit_a_err = fit.ParError(0)
    fit_b = fit.Value(1)
    fit_b_err = fit.ParError(1)

    print("fit: ",fit_a, fit_b)

    # line with statistical error
    xx_fit = np.linspace(xmin,xmax,50)
    yy_fit = fit_a + xx_fit * fit_b
    yy_fit_err_up = fit_a+fit_a_err + xx_fit * (fit_b+fit_b_err)
    yy_fit_err_down = fit_a-fit_a_err + xx_fit * (fit_b-fit_b_err)
    fit_line = ROOT.TGraphAsymmErrors(len(xx_fit), xx_fit, yy_fit, np.zeros(len(xx_fit)), np.zeros(len(xx_fit)), yy_fit - yy_fit_err_down, yy_fit_err_up - yy_fit)
    fit_line.SetFillColorAlpha(ROOT.kAzure-4,0.4)
    fit_line.SetLineColor(ROOT.kBlack)

    canvas=ROOT.TCanvas("canvas","canvas",1000,600)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.03)
    canvas.SetTopMargin(0.015)
    canvas.SetBottomMargin(0.15)
    canvas.SetGrid()

    textsize = 24./(canvas.GetWh()*canvas.GetAbsHNDC())

    entry = ROOT.TLegendEntry()
    entry.SetLabel("Fit (#pm stat.)")
    entry.SetOption("fl")
    entry.SetTextSize(textsize)
    entry.SetTextFont(42)
    entry.SetLineColor(ROOT.kBlack)
    #entry.SetLineWidth(2)
    entry.SetFillColorAlpha(ROOT.kAzure-4,0.4)
    entry.SetTextAlign(12)
    entry.SetMarkerSize(0)

    legend=ROOT.TLegend(0.2,0.8,0.5,0.9)
    legend.AddEntry(graphXsecL, "2017 High PU (#pm stat.)", "pe")
    if plot_lowPU:
        legend.GetListOfPrimitives().Add(entry2)

    legend=ROOT.TLegend(0.16,0.85,0.4,0.97)
    legend.AddEntry(graphXsecL, "Measurements (#pm stat.)", "pe")
    legend.GetListOfPrimitives().Add(entry)
    legend.SetTextSize(textsize)
    legend.SetTextFont(42)

    graphXsecL.GetYaxis().SetTitleFont(42)
    graphXsecL.GetYaxis().SetTitleSize(textsize*1.2)
    graphXsecL.GetYaxis().SetLabelFont(42)
    graphXsecL.GetYaxis().SetLabelSize(textsize*1.2)
    graphXsecL.GetXaxis().SetTitleFont(42)
    graphXsecL.GetXaxis().SetTitleSize(textsize*1.2)
    graphXsecL.GetXaxis().SetLabelFont(42)
    graphXsecL.GetXaxis().SetLabelSize(textsize*1.2)

    graphXsecL.Draw("AP")

    fit_line.Draw("3L same")
    legend.Draw("same")
    # latex.SetTextAlign(31)

    # latex.DrawLatex(0.97, 0.95, "2017")
    latex.SetTextFont(42)

    if region == 'BB':
        str = "barrel-barrel"
    elif region == 'BE':
        str = "barrel-endcap"
    elif region == 'EE':
        str = "endcap-endcap"
    else:
        str = "inclusive"

    latex.SetTextSize(textsize)
    latex.DrawLatex(0.77, 0.93, year)
    latex.DrawLatex(0.77, 0.875, str)
    # latex.DrawLatex(0.12, 0.95, title)

    latex.SetTextSize(textsize*1.2)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.22, 0.18, "Preliminary")

    latex.SetTextAlign(11)
    latex.SetTextFont(62)
    latex.DrawLatex(0.15, 0.18, 'CMS')

    outstring = "{0}/{1}_vs_{2}".format(outDir,ZXSec,xAxis)
    if run_range:
        outstring += "_run{0}to{1}".format(*run_range)

    canvas.SaveAs(outstring+".png")
    canvas.Close()

# ------------------------------------------------------------------------------
print(">>> load csv file in dataframe")
data = pd.read_csv(args.rates, sep=',')[
    ['fill','run','lumiRec','timewindow','pileUp','tdate_begin','tdate_end',
    'zBB_relstat', 'zBE_relstat', 'zEE_relstat',
    'xsecBB','xsecBE','xsecEE',
    'xsecBB_mc','xsecBE_mc','xsecEE_mc',
    ]]

data['xsec'] = data['xsecBB'] + data['xsecBE'] + data['xsecEE']
data['xsec_mc'] = data['xsecBB_mc'] + data['xsecBE_mc'] + data['xsecEE_mc']

data['z_relstat'] = np.sqrt(data['zEE_relstat']**2 + data['zBB_relstat']**2) + data['zBE_relstat']

# make_plots(data, "xsec",      'z_relstat',   title="Z rate / lumi - 2017 All")
#make_plots(data, "xsec_mc",   'z_relstat',   title="Z rate / lumi - 2017 All", normalized=True)
# make_plots(data, "xsecBB_mc", 'zBB_relstat', title="corrected Z BB rate / lumi - 2017 All", normalized=True)
# make_plots(data, "xsecBE_mc", 'zBE_relstat', title="corrected Z BE rate / lumi - 2017 All", normalized=True)
# make_plots(data, "xsecEE_mc", 'zEE_relstat', title="corrected Z EE rate / lumi - 2017 All", normalized=True)
#

### 2017

# make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', year="2017 B", title="corrected", run_range=(297046,299329), normalized=True)
# make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', year="2017 B", title="corrected", run_range=(297046,299329), normalized=True)
# make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', year="2017 B", title="corrected", run_range=(297046,299329), normalized=True)
#
# make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', year="2017 C", title="corrected", run_range=(299368,302029), normalized=True)
# make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', year="2017 C", title="corrected", run_range=(299368,302029), normalized=True)
# make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', year="2017 C", title="corrected", run_range=(299368,302029), normalized=True)
#
# make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', year="2017 D", title="corrected", run_range=(302030,303434), normalized=True)
# make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', year="2017 D", title="corrected", run_range=(302030,303434), normalized=True)
# make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', year="2017 D", title="corrected", run_range=(302030,303434), normalized=True)
#
# make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', year="2017 E", title="corrected", run_range=(303434,304797), normalized=True)
# make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', year="2017 E", title="corrected", run_range=(303434,304797), normalized=True)
# make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', year="2017 E", title="corrected", run_range=(303434,304797), normalized=True)
#
# make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', year="2017 F", title="corrected", run_range=(305040,306462), normalized=True)
# make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', year="2017 F", title="corrected", run_range=(305040,306462), normalized=True)
# make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', year="2017 F", title="corrected", run_range=(305040,306462), normalized=True)
#
# ### total 2017
# make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', year="2017", title="corrected", run_range=(297046,306462), normalized=True)
# make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', year="2017", title="corrected", run_range=(297046,306462), normalized=True)
# make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', year="2017", title="corrected", run_range=(297046,306462), normalized=True)
#

### 2018
make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', title="corrected", year="2018 A", run_range=(315252,316995), normalized=True)
make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', title="corrected", year="2018 A", run_range=(315252,316995), normalized=True)
make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', title="corrected", year="2018 A", run_range=(315252,316995), normalized=True)

make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', title="corrected", year="2018 B", run_range=(317080,319310), normalized=True)
make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', title="corrected", year="2018 B", run_range=(317080,319310), normalized=True)
make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', title="corrected", year="2018 B", run_range=(317080,319310), normalized=True)

make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', title="corrected", year="2018 C", run_range=(319337,320065), normalized=True)
make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', title="corrected", year="2018 C", run_range=(319337,320065), normalized=True)
make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', title="corrected", year="2018 C", run_range=(319337,320065), normalized=True)

## total 2018
make_plots(data, "xsecBB_mc", 'zBB_relstat', region='BB', title="corrected", year="2018", run_range=(315252,320065), normalized=True)
make_plots(data, "xsecBE_mc", 'zBE_relstat', region='BE', title="corrected", year="2018", run_range=(315252,320065), normalized=True)
make_plots(data, "xsecEE_mc", 'zEE_relstat', region='EE', title="corrected", year="2018", run_range=(315252,320065), normalized=True)

# make_plots(data, "xsec_mc",   title="Corrected Z rate / lumi - 2017 All", plot_fills=True)
# make_plots(data, "zXSec_mc", title="Corrected Z rate / lumi - 2017B-E", run_range=(297046,304797))
# make_plots(data, "zXSec_mc", title="Corrected Z rate / lumi - 2017B", run_range=(297046,299329))
# make_plots(data, "xsec",      title="Z rate / lumi - 2017C", run_range=(299368,302029))
# make_plots(data, "xsec",      title="Z rate / lumi - 2017D", run_range=(302030,303434))
# make_plots(data, "xsec",      title="Z rate / lumi - 2017E", run_range=(303434,304797))
# make_plots(data, "xsec_mc",   title="Corrected Z rate / lumi - 2017C", run_range=(299368,302029))
# make_plots(data, "xsec_mc",   title="Corrected Z rate / lumi - 2017D", run_range=(302030,303434))
# make_plots(data, "xsec_mc",   title="Corrected Z rate / lumi - 2017E", run_range=(303434,304797))
# make_plots(data, "zXSec_mc", title="Corrected Z rate / lumi - 2017F", run_range=(305040,306462))

# make_plots(data, "zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi", nround=4)
# make_plots(data, "zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017B", run_range=(297046,299329), nround=4)
# make_plots(data, "zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017C", run_range=(299368,302029), nround=4)
# make_plots(data, "zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017D", run_range=(302030,303434), nround=4)
# make_plots(data, "zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017E", run_range=(303434,304797), nround=4)
# make_plots(data, "zXSec_mc", xAxis='time', xTitle="time", title="Corrected Z rate / lumi - 2017F", run_range=(305040,306462), nround=4)
