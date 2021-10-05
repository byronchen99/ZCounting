import ROOT
import pandas as pd
import numpy as np
import argparse
import pdb
import os, sys
import math
import uncertainties as unc
from scipy.optimize import curve_fit
from scipy import stats

sys.path.append(os.getcwd())
print(os.getcwd())

from python.corrections import apply_muon_prefire, apply_ECAL_prefire, apply_pileup_correction

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime, cms, preliminary, text, workinprogress, unorm, linear

pd.options.mode.chained_assignment = None

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)


latex = ROOT.TLatex()
latex.SetNDC()

parser = argparse.ArgumentParser()

parser.add_argument("-r","--rates", required=True, type=str, help="csv file with z rates per Measurement")
parser.add_argument("-x","--xsec",  type=str, help="csv file where xsec should be taken from (e.g. from low pileup run)")
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
    xsecBB = sum(data_xsec['zDelBB'])/sum(data_xsec['lumiRec'])
    xsecBE = sum(data_xsec['zDelBE'])/sum(data_xsec['lumiRec'])
    xsecEE = sum(data_xsec['zDelEE'])/sum(data_xsec['lumiRec'])

textsize = 0.03

# ------------------------------------------------------------------------------
def make_plots(df,
    yAxis,
    yLabel="sigma",
    xAxis='lumi',
    run_range=None,
    title="",
    year="2017",
    region='inclusive',
    plot_fills=False,
    plot_all=False,
    nround=0,
    normalized=False,
    resource=""
):
    """
    valid xAxis: 'lumi', 'pileUp', 'measurement', 'time'
    """

    if run_range:
        data = df.loc[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])]
        if len(df) ==0:
            return
    else:
        data = df

    data['y_Err'] = data[yAxis].apply(lambda x: x.s)
    data['y'] = data[yAxis].apply(lambda x: x.n)
    if yAxis.replace("_mc","") in data.keys():
        data['y0_Err'] = data[yAxis.replace("_mc","")].apply(lambda x: x.s)
        data['y0'] = data[yAxis.replace("_mc","")].apply(lambda x: x.n)
    else:
        data['y0_Err'] = data['y_Err']
        data['y0'] = data['y']

    if sum(data['y'].isnull()) > 0:
        print(">>> sort out {0} points with nan".format(sum(data['y'].isnull())))
        data = data.loc[~data['y'].isnull()]

    data['relstat'] = data['y_Err'] / data['y']

    # x_step: intervall in which the measurements are collected in one point

    if xAxis == 'lumi':
        xAxis = 'lumiRecInst'
        data[xAxis] = data['lumiRec'] / data['timewindow']
        xTitle="instantaneous recorded luminosity [Hz/pb]"
        x_step = 0.0005
    elif xAxis == 'time':
        data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2
        #data['time'] -= 7200
        data['time'] -= (ROOT.TDatime(2042,01,01,12,00,00).Convert() - ROOT.TDatime(2017,01,01,12,00,00).Convert())
        xTitle="time"
        x_step = 1
    elif xAxis == 'pileUp':
        xTitle="average number of pileup interactions"
        x_step = 2
    else:
        xTitle=xAxis
        x_step = 0.0005

    if yLabel == "sigma":
        yLabel = "#sigma^{fid}_{Z} [pb]"

        if normalized:
            data[yAxis] /= np.mean(data[yAxis].values)
            yLabel = "#sigma^{fid}_{Z} / <#sigma^{fid}_{Z}>"

    elif resource != "":
        # plot an efficiency on y axis - load necessery resources

        # >>>  load MC info and format
        if len(year.split(" ")) == 1:
            suffix = year
        elif len(year.split(" ")) == 0:
            print("Invalid year")
        elif year.split(" ")[0] == "2016":
            if year == "2016 pre VFP" or year.split(" ")[1] in ("B","C","D","E","F"):
                suffix="2016preVFP"
            elif year == "2016 post VFP" or year.split(" ")[1] in ("G","H"):
                suffix="2016postVFP"
        else:
            suffix = year.split(" ")[0]

        genInfoFile = 'res/GenInfo-V11_15_03-d20210611-t120445/infoMC_gen_{0}.json'.format(suffix)
        recoInfoFile = 'res/RecoInfo-V11_15_03-d20210611-t120433/infoMC_reco_{0}.json'.format(suffix)

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

    data['y_Err'] = data[yAxis].apply(lambda x: x.s)
    data['y'] = data[yAxis].apply(lambda x: x.n)
    if yAxis.replace("_mc","") in data.keys():
        data['y0_Err'] = data[yAxis.replace("_mc","")].apply(lambda x: x.s)
        data['y0'] = data[yAxis.replace("_mc","")].apply(lambda x: x.n)
    else:
        data['y0_Err'] = data['y_Err']
        data['y0'] = data['y']

    data['relstat'] = data['y_Err'] / data['y']


    # print(">>> sort out {0} points with low statistic".format(len(data['relstat'] > 0.05)))
    # data = data.loc[data["relstat"] < 0.05]

    if plot_fills:
        print(">>> make plots for each fill")
        # For each Fill
        for fill in data.drop_duplicates('fill')['fill'].values:
            dFill = data.loc[data['fill'] == fill]

            subDir = outDir+"/PlotsXSec_"+str(fill)
            if not os.path.isdir(subDir):
                os.mkdir(subDir)

            graphXsecL=ROOT.TGraphErrors(len(dFill),dFill[xAxis].values,dFill['y'].values,np.zeros(len(dFill)),data['y_Err'].values)
            graphXsecL.SetName("graph_metaXsecAtlas")
            graphXsecL.SetMarkerStyle(23)
            graphXsecL.SetMarkerColor(ROOT.kAzure-4)
            graphXsecL.SetTitle(title)
            graphXsecL.GetXaxis().SetTitle(xTitle)
            graphXsecL.GetYaxis().SetTitle(yLabel)

            graphXsecL.GetYaxis().SetTitleOffset(1.0)

            print(">>> cross sections for Fill "+str(fill))
            print(">>> the simple average cross section is "+str(sum(dFill['y'].values)/len(dFill)))
            graphXsecL.Fit("pol1","","",0.0,0.025)
            c3=ROOT.TCanvas("c3","c3",1000,600)
            c3.SetGrid()
            graphXsecL.Draw("AP")


            legend=ROOT.TLegend(0.2,0.8,0.4,0.9)
            legend.AddEntry(graphXsecL,"Measurement (#pm stat.)","pe")
            legend.Draw("same")

            c3.SaveAs(subDir+"/"+yAxis+"_vs_"+xAxis+".png")
            c3.Delete()

    if plot_all:
        print(">>> make plot all measurements")
        if xAxis == 'measurement':
            xPoints = np.arange(0,len(data),1.)
        else:
            xPoints = data[xAxis].values

        graphXsecL=ROOT.TGraphErrors(len(data), xPoints, data['y'].values, np.zeros(len(data)), data['y_Err'].values)
        graphXsecL.SetName("graph_metaXsecAtlas")
        graphXsecL.SetMarkerStyle(33)
        graphXsecL.SetMarkerColor(ROOT.kAzure-4)
        graphXsecL.SetMarkerSize(1.)
        graphXsecL.SetTitle(title)
        graphXsecL.GetXaxis().SetTitle(xTitle)
        graphXsecL.GetYaxis().SetTitle(yTitle)
        graphXsecL.GetYaxis().SetTitleOffset(1.0)

        #graphXsecL.GetYaxis().SetRangeUser(500.,700.)

        print(">>> Producing all cross sections")
        print(">>> the simple average cross section is "+str(sum(data['y'].values)/len(data)))
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

        fit = graphXsecL.Fit("pol1","0S","",xmin, xmax)
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

        outstring = "{0}/{1}_vs_{2}".format(outDir,yAxis,xAxis)
        if run_range:
            outstring += "_run{0}to{1}".format(*run_range)

        c3.SaveAs(outstring+".png")
        c3.Close()

    print(">>> make plot and combine measurements into bins")
    xx = np.arange(data[xAxis].min()*(1./x_step)//1/(1./x_step),data[xAxis].max()*(1./x_step)//1/(1./x_step), x_step)
    xx = np.append(xx, xx[-1]+x_step)
    xx_centers = []
    yy = []
    yy_err = []
    yy0 = []
    yy0_err = []
    for i in range(0,len(xx)-1):
        dyy = data.loc[(data[xAxis] < xx[i+1]) & (data[xAxis] >= xx[i])][[xAxis, 'y', 'y_Err', 'y0', 'y0_Err']]

        if len(dyy)==0:
            continue

        yy_avg = np.average(dyy['y'].values, weights=1./dyy['y_Err'].values)
        yy_avg_err = np.sqrt(1./sum((1./dyy['y_Err'].values)**2))

        yy0_avg = np.average(dyy['y0'].values, weights=1./dyy['y0_Err'].values)
        yy0_avg_err = np.sqrt(1./sum((1./dyy['y0_Err'].values)**2))

        if yy_avg_err/yy_avg > 0.05:
            continue

        xx_centers.append((xx[i] + (xx[i+1] - xx[i]) / 2.))
        yy.append(yy_avg)
        yy_err.append(yy_avg_err)
        yy0.append(yy0_avg)
        yy0_err.append(yy0_avg_err)

    graphXsecL=ROOT.TGraphErrors(len(xx_centers), np.array(xx_centers), np.array(yy), np.ones(len(xx_centers))*x_step/2, np.array(yy_err))
    graphXsecL.SetName("graph_metaXsecAtlas")
    graphXsecL.SetMarkerStyle(33)
    graphXsecL.SetMarkerColor(ROOT.kAzure-4)
    graphXsecL.SetMarkerSize(1.)
    graphXsecL.SetTitle("")
    graphXsecL.GetXaxis().SetTitle(xTitle)
    graphXsecL.GetYaxis().SetTitle(yLabel)
    graphXsecL.GetYaxis().SetTitleOffset(1.2)

    if '_mc' in yAxis:
        graphXsecL0=ROOT.TGraphErrors(len(xx_centers), np.array(xx_centers), np.array(yy0), np.ones(len(xx_centers))*x_step/2, np.array(yy0_err))
        graphXsecL0.SetName("graph_metaXsecAtlas")
        graphXsecL0.SetMarkerStyle(33)
        graphXsecL0.SetMarkerColor(ROOT.kGray)
        graphXsecL0.SetLineColor(ROOT.kGray)
        graphXsecL0.SetMarkerSize(1.)
        graphXsecL0.SetTitle("")
        graphXsecL0.GetXaxis().SetTitle(xTitle)
        graphXsecL0.GetYaxis().SetTitle(yLabel)
        graphXsecL0.GetYaxis().SetTitleOffset(1.2)

    xmin = min(xx_centers) - x_step
    xmax = max(xx_centers) + x_step
    ymin = min(np.array(yy)-np.array(yy_err))
    xRange = (xmin, xmax)

    if "sigma" in yLabel:
        graphXsecL.SetMarkerColor(ROOT.kAzure-4)
        ymax = max(np.array(yy)+np.array(yy_err))
        yRange = (ymin - 0.15*(ymax-ymin), ymax + 0.25*(ymax-ymin))
    elif resource != "":
        yTrue = iTrue[resource].values
        yTrueErr = iTrue[resource+"_err"].values
        xTrue = iTrue['center'].values

        yReco = iReco[resource].values
        yRecoErr = iReco[resource+"_err"].values
        xReco = iReco['center'].values

        gTrue=ROOT.TGraphErrors(len(yTrue), xTrue, yTrue, np.ones(len(xTrue))*x_step/2, np.array(yRecoErr))
        gTrue.SetLineColor(ROOT.kRed)
        # gTrue.SetLineWidth(2)
        # gTrue.SetLineStyle(1)
        gTrue.SetMarkerColor(2)
        gTrue.SetMarkerStyle(20)
        gTrue.SetMarkerSize(0.5)
        gTrue.SetTitle("")

        gReco=ROOT.TGraphErrors(len(xReco), xReco, yReco, np.ones(len(xReco))*x_step/2, np.array(yRecoErr))

        gReco.SetLineColor(ROOT.kOrange)
        # gReco.SetLineWidth(2)
        # gReco.SetLineStyle(2)
        gReco.SetMarkerColor(ROOT.kOrange)
        gReco.SetMarkerStyle(21)
        gReco.SetMarkerSize(0.5)
        gReco.SetTitle("")

        yMax = max(max(yTrue), max(yReco))
        ymax = max(max(np.array(yy)+np.array(yy_err)),yMax)
        yRange = (ymin - 0.15*(ymax-ymin), ymax + 0.35*(ymax-ymin))

        graphXsecL.SetMarkerColor(1)


    graphXsecL.GetXaxis().SetRangeUser(*xRange)
    graphXsecL.GetYaxis().SetRangeUser(*yRange)

    def fit(func):
        popt, pcov = curve_fit(func, np.array(xx_centers), np.array(yy), sigma=np.array(yy_err))
        perr = np.sqrt(np.diag(pcov))

        params = unc.correlated_values(popt, pcov)
        yyFit = func(np.array(xx_centers), *params)
        yyFit_err = np.array([iy.s for iy in yyFit])
        yyFit_val = np.array([iy.n for iy in yyFit])

        chi2 = sum(((yy - yyFit_val) / yyFit_err)**2)
        ndf = len(xx_centers) - len(popt)

        p_chi2 = stats.distributions.chi2.sf(chi2, ndf)
        return p_chi2, params, func

    lChi2, params, func = fit(linear)
    funcname = "Linear fit"
    # qChi2, qParams, qFunc = fit(quad)
    # if lChi2 > qChi2:
    #     funcname = "Linear fit"
    #     params = lParams
    #     func = lFunc
    # else:
    #     funcname = "Quadratic fit"
    #     params = qParams
    #     func = qFunc

    # line with statistical error
    xx_fit = np.linspace(xmin,xmax,50)
    yy_fit = np.array([x.n for x in func(xx_fit, *params)])
    yy_fit_err_up   = np.array([x.s for x in func(xx_fit, *params)])
    yy_fit_err_down = np.array([x.s for x in func(xx_fit, *params)])
    fit_line = ROOT.TGraphAsymmErrors(len(xx_fit), xx_fit, yy_fit, np.zeros(len(xx_fit)), np.zeros(len(xx_fit)), yy_fit_err_down, yy_fit_err_up)
    fit_line.SetFillColorAlpha(ROOT.kAzure-4,0.4)
    fit_line.SetLineColor(ROOT.kBlack)

    canvas=ROOT.TCanvas("canvas","canvas",1000,600)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.03)
    canvas.SetTopMargin(0.015)
    canvas.SetBottomMargin(0.15)
    canvas.SetGrid()

    textsize = 24./(canvas.GetWh()*canvas.GetAbsHNDC())

    graphXsecL.GetYaxis().SetTitleFont(42)
    graphXsecL.GetYaxis().SetTitleSize(textsize*1.2)
    graphXsecL.GetYaxis().SetLabelFont(42)
    graphXsecL.GetYaxis().SetLabelSize(textsize*1.2)
    graphXsecL.GetXaxis().SetTitleFont(42)
    graphXsecL.GetXaxis().SetTitleSize(textsize*1.2)
    graphXsecL.GetXaxis().SetLabelFont(42)
    graphXsecL.GetXaxis().SetLabelSize(textsize*1.2)

    entry = ROOT.TLegendEntry()
    entry.SetLabel("{0} (#pm stat.)".format(funcname))
    entry.SetOption("fl")
    entry.SetTextSize(textsize)
    entry.SetTextFont(42)
    entry.SetLineColor(ROOT.kBlack)
    #entry.SetLineWidth(2)
    entry.SetFillColorAlpha(ROOT.kAzure-4,0.4)
    entry.SetTextAlign(12)
    entry.SetMarkerSize(0)

    if "sigma" in yLabel:
        legend=ROOT.TLegend(0.72,0.82,0.97,0.97)
        legend.GetListOfPrimitives().Add(entry)
    else:
        # legend=ROOT.TLegend(0.16,0.80,0.4,0.97)
        # legend.AddEntry(gReco, "MC eff (T&P)", "pe")
        legend=ROOT.TLegend(0.72,0.82,0.97,0.97)
        legend.AddEntry(gTrue, "MC", "pe")

    graphXsecL.Draw("AP")

    if '_mc' in yAxis:
        graphXsecL0.Draw("pe same")
        legend.AddEntry(graphXsecL0, "Measurements w/o correction", "pe")

    legend.AddEntry(graphXsecL, "Measurements", "pe")
    legend.SetTextSize(textsize)
    legend.SetTextFont(42)


    if "sigma" in yLabel:
        fit_line.Draw("3L same")
    else:
        # gReco.Draw("pe same")
        gTrue.Draw("pe same")

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
    elif region == "B":
        str = "barrel"
    elif region == "E":
        str = "endcap"
    else:
        str = "inclusive"

    latex.SetTextSize(textsize)
    latex.DrawLatex(0.15, 0.23, year)
    latex.DrawLatex(0.15, 0.18, str)
    # latex.DrawLatex(0.12, 0.95, title)

    cms(x=0.15, y=0.93, textsize=textsize)
    workinprogress(x=0.23, y=0.93, textsize=textsize)

    outstring = "{0}/{1}_vs_{2}".format(outDir,yAxis,xAxis)
    outstring += year.replace(" ","_")

    # if run_range:
    #     outstring += "_run{0}to{1}".format(*run_range)

    canvas.SaveAs(outstring+".png")
    canvas.SaveAs(outstring+".pdf")
    canvas.Close()

# ------------------------------------------------------------------------------
print(">>> load csv file in dataframe")

data_rates = pd.read_csv(args.rates, sep=',')[
    ['fill','run','lumiRec','timewindow','pileUp','tdate_begin','tdate_end',
    'zDelBB_mc','zDelBE_mc','zDelEE_mc',
    'zDelBB','zDelBE','zDelEE',
    'ZBBeff','ZBEeff','ZEEeff',
    'HLTeffBB', 'HLTeffBE', 'HLTeffEE',
    'SeleffB', 'SeleffE',
    'GloeffB', 'GloeffE',
    'TrkeffB', 'TrkeffE',
    'StaeffB', 'StaeffE',
    # 'ZBBeff_mc','ZBEeff_mc','ZEEeff_mc',
    ]]

# convert to uncertainties
data_rates['HLTeffBB'] = data_rates['HLTeffBB'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['HLTeffBE'] = data_rates['HLTeffBE'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['HLTeffEE'] = data_rates['HLTeffEE'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['SeleffB'] = data_rates['SeleffB'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['SeleffE'] = data_rates['SeleffE'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['TrkeffB'] = data_rates['TrkeffB'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['TrkeffE'] = data_rates['TrkeffE'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['StaeffB'] = data_rates['StaeffB'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['StaeffE'] = data_rates['StaeffE'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['GloeffB'] = data_rates['GloeffB'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['GloeffE'] = data_rates['GloeffE'].apply(lambda x: unc.ufloat_fromstr(x))

data_rates['ZBBeff'] = data_rates['ZBBeff'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['ZBEeff'] = data_rates['ZBEeff'].apply(lambda x: unc.ufloat_fromstr(x))
data_rates['ZEEeff'] = data_rates['ZEEeff'].apply(lambda x: unc.ufloat_fromstr(x))

# take uncertainties from uncorrected zDel
data_rates['zDelBB_mc'] = data_rates['zDelBB'].apply(lambda x: unc.ufloat_fromstr(x)/unc.ufloat_fromstr(x).n) * data_rates['zDelBB_mc']
data_rates['zDelBE_mc'] = data_rates['zDelBE'].apply(lambda x: unc.ufloat_fromstr(x)/unc.ufloat_fromstr(x).n) * data_rates['zDelBE_mc']
data_rates['zDelEE_mc'] = data_rates['zDelEE'].apply(lambda x: unc.ufloat_fromstr(x)/unc.ufloat_fromstr(x).n) * data_rates['zDelEE_mc']

apply_muon_prefire(data_rates)
apply_ECAL_prefire(data_rates)
# apply_pileup_correction(data_rates)

data_rates['xsecBB'] = data_rates['zDelBB_mc'] / data_rates['lumiRec']
data_rates['xsecBE'] = data_rates['zDelBE_mc'] / data_rates['lumiRec']
data_rates['xsecEE'] = data_rates['zDelEE_mc'] / data_rates['lumiRec']
data_rates['xsec'] = data_rates['xsecBB'] + data_rates['xsecBE'] + data_rates['xsecEE']

for yy, ylabel, region, mcRes in (
    ("xsec", "sigma", "", ""),
    # ("xsecBB", "sigma", "BB", ""),
    # ("xsecBE", "sigma", "BE", ""),
    # ("xsecEE", "sigma", "EE", ""),
    # ("xsec_mc", "sigma", "", ""),
    # ("xsecBB_mc", "sigma", "BB", ""),
    # ("xsecBE_mc", "sigma", "BE", ""),
    # ("xsecEE_mc", "sigma", "EE", ""),
    # ("ZBBeff", "Z efficiency","BB", "effBB"),
    # ("ZBEeff", "Z efficiency","BE", "effBE"),
    # ("ZEEeff", "Z efficiency","EE", "effEE"),
    # ('HLTeffBB' ,'Muon HLT efficiency', "BB", "HLTB"),
    # ('HLTeffBE' ,'Muon HLT efficiency', "BE", "HLTE"),
    # ('HLTeffEE' ,'Muon HLT efficiency', "EE", "HLTE"),
    # ('SeleffB' ,'Muon selelction efficiency', "B", "SelB"),
    # ('SeleffE' ,'Muon selelction efficiency', "E", "SelE"),
    # ('GloeffB' ,'Global muon efficiency', "B", "GloB"),
    # ('GloeffE' ,'Global muon efficiency', "E", "GloE"),
    # ('TrkeffB' ,'Muon inner track efficiency', "B", "TrkB"),
    # ('TrkeffE' ,'Muon inner track efficiency', "E", "TrkE"),
    # ('StaeffB' ,'Muon standalone efficiency', "B", "StaB"),
    # ('StaeffE' ,'Muon standalone efficiency', "E", "StaE"),
):
    # # single eras
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 B", run_range=(272007,275376), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 C", run_range=(275657,276283), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 D", run_range=(276315,276811), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 E", run_range=(276831,277420), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 F", run_range=(277772,278808), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 G", run_range=(278820,280385), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 H", run_range=(280919,284044), normalized=False)

    # 2017
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 B", title="corrected", run_range=(297046,299329), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 C", title="corrected", run_range=(299368,302029), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 D", title="corrected", run_range=(302030,303434), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 E", title="corrected", run_range=(303434,304797), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017 F", title="corrected", run_range=(305040,306462), normalized=False)

    # ## 2018
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 A", run_range=(315252,316995), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 B", run_range=(317080,319310), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 C", run_range=(319337,320065), normalized=False)
    make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 D", run_range=(320673,325175), normalized=False)

    # # total 2016 pre VFP
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 pre VFP", run_range=(272007,278769), normalized=True)
    # # total 2016 post VFP
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 post VFP", run_range=(278769,294645), normalized=True)
    # ## total 2017
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, year="2017", title="corrected", run_range=(297046,306462), normalized=True)
    # ## total 2018
    # make_plots(data_rates, yAxis=yy, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018", run_range=(315252,325175), normalized=True)
