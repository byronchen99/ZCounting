###
#   Script to illustrate background functions used for fitting
#
###

import ROOT
import os
import numpy as np
import json
import pdb
import uncertainties as unc

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))

# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

### Settings ###
# ---------------------->>>
massLo = 60
massHi = 120
massBin = 120
ptCut = 25
etaCut = 2.4
# <<<------------------------


################################################################################
if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", help="Output directory to the histogram file",
                        default="./Test")
    parser.add_argument("-l", "--logscale", help="Plot in logscale",
                        default=False, action="store_true")

    args = parser.parse_args()
    logscale = args.logscale
    output = args.output

    print("INFO: Loading C marco...")
    # load functions for fitting
    ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(
        __file__)) + "/calculateDataEfficiency.C")

    ROOT.set_massRange(massLo-5., massHi+5., massBin)
    ROOT.set_ptCut(ptCut)
    ROOT.set_etaCut(etaCut)
    ROOT.set_output(output)    

    # # fill histogram with gaussian distribution
    # ROOT.TH1D hPV = ROOT.TH1F("hPV", "hPV", 200, -5,5);
    # TF1 *f1 = new TF1("f1", "[2]*TMath::Gaus(x,[0],[1])");
    # f1->SetParameters(1,1,1);
    # h1->FillRandom("f1");
    # h1->Draw();
    print("Generate template")
    templates = ROOT.generateTemplate_ZYield("/afs/desy.de/user/d/dwalter/nfsHome/data/Lumi/V17_02/TTrees/ZCountingAll-V17_02-Fall17-DYJetsToLL_M_50_LO.root",0,0)
    # print("load histogram templates")
    # templates = ROOT.TFile("histTemplates_HLT.root","READ")
    
    h2 = templates.Get("h_mass_2hlt_I")
    h1 = templates.Get("h_mass_2hlt_I")

    h2.SetDirectory(0)
    h1.SetDirectory(0)
    templates.Close()
    templates.Delete()

    if not os.path.isdir(output):
        os.mkdir(output)

    canvas=ROOT.TCanvas("canvas","canvas",1000,600)
    margin_left = 0.125
    margin_right = 0.03
    margin_top = 0.015
    margin_bottom = 0.15

    canvas.SetLeftMargin(margin_left)
    canvas.SetRightMargin(margin_right)
    canvas.SetTopMargin(margin_top)
    canvas.SetBottomMargin(margin_bottom)
    canvas.SetGrid()

    textsize = 24./(canvas.GetWh()*canvas.GetAbsHNDC())

    param_mass = ROOT.RooRealVar("m", "mass", massLo-5., massHi+5.)

    mframe = ROOT.RooPlot("mframe","m frame", param_mass, massLo, massHi, massBin)
    if logscale:
        # legend = ROOT.TLegend(margin_left+0.025, margin_bottom+0.025, 1-margin_right-0.05, margin_bottom+0.35)
        legend = ROOT.TLegend(margin_left+0.025, margin_bottom+0.025, 1-margin_right-0.5, margin_bottom+0.25)
        legend.SetNColumns(2)
    else:
        legend = ROOT.TLegend(margin_left+0.025, 0.5, 0.45, 1-margin_top-0.01)
        legend.SetNColumns(1)

    legend.SetBorderSize(0)

    functions = []
    entries=[]
    for i, args in enumerate([
    #    ("BWxCB", 0.2, 1.1, 1.0, 2.6),
       # ("BWxCB", 0.1, 1.2, 1.3, 4.2),
       # ("BWxCB", -0.6, 3, 18, 0.5) ,
       # ("BWxCB", 0.8, 0.1, 0.1, 10) ,
       # ("BWxCB", 0.2, 0.9, 1.5, 2.3) ,
       # ("BWxCB", 0.1, 1.2, 1.2, 3.7) ,
        ("MCxGauss", 0, 1),
        ("MCxGauss", 0, 4),
        ("MCxGauss", 0, 0.2),
        ("MCxGauss", 2, 1),
       # ("MC", 0, ) ,
       # ("BW", 0, ) ,
    #    ("BWxCB", 0, 1, 5, 5) ,
       
    ]):
        model = args[0]
        params = args[1:]

        entry = ROOT.TLegendEntry()
        entry.SetTextSize(textsize)
        
        if model == "BWxCB":
            fX = ROOT.CBreitWignerConvCrystalBall(param_mass, 0, 0)
            fX.mean.setVal(params[0])
            fX.sigma.setVal(params[1])
            fX.alpha.setVal(params[2])
            fX.n.setVal(params[3])

            entry.SetLabel("#mu={0}, #sigma={1}".format(params[0],params[1])
            +", #alpha="+str(params[2])
            +", n="+str(params[3]))
                    
        elif model == "MCxGauss":
            fX = ROOT.CMCTemplateConvGaussian(param_mass, h2, 0, 0)
            fX.mean.setVal(params[0])
            fX.sigma.setVal(params[1])
            entry.SetLabel("#mu={0}, #sigma={1}".format(params[0],params[1]))
        elif model == "MC":
            #inHist = ROOT.TH1D(h2.Clone("h2MC"))
            #dataHist = ROOT.RooDataHist("h2MC","h2MC",ROOT.RooArgSet(param_mass),inHist)
            fX  = h2#ROOT.RooHistPdf("h2MC","h2MC",param_mass,dataHist,0)
            entry.SetLabel("MC template".format(params[0],params[1]))
        elif model == "BW":
            fX = ROOT.RooBreitWigner("BW","BW",param_mass,91.1876,2.4952)
            entry.SetLabel("BreitWigner".format(params[0],params[1]))
            
        entry.SetOption("L")
        entry.SetLineStyle(1)
        entry.SetLineColor(i+1)
        entry.SetLineWidth(2)
        entry.SetTextAlign(12)
        legend.GetListOfPrimitives().Add(entry)
        entries.append(entry)

        # fX.t1.setVal(params[4])
        # fX.frac.setVal(params[5])
        fX.model.plotOn(mframe)
        mframe.getAttLine().SetLineColor(i+1)
        functions.append(fX)

    if logscale:
        canvas.SetLogy(1)
        mframe.SetMinimum(0.000001)
        mframe.SetMaximum(mframe.GetMaximum()*2);

    mframe.SetTitle("")
    mframe.SetYTitle("a.u.")
    mframe.SetXTitle("mass [GeV]")
    mframe.SetTitleOffset(1.1, "Y")
    # mframe.SetMaximum(0.06)
    mframe.Draw()

    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(11)
    latex.SetTextFont(42)
    latex.SetTextSize(textsize*1.5)
    # latex.DrawLatex(0.2, 0.9, model)
    latex.DrawLatex(0.7, 0.85, model)
    
    saveas = output+"/"+model
    if logscale:
        saveas += "_logscale"
    canvas.SaveAs(saveas+".png")
    canvas.SaveAs(saveas+".eps")
    canvas.Close()
