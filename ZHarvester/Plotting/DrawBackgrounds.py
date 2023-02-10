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
# <<<------------------------


################################################################################
if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", help="Output directory to the histogram file",
                        default="./Test")
    args = parser.parse_args()

    output = args.output

    print("INFO: Loading C marco...")
    # load functions for fitting
    ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(
        __file__)) + "/calculateDataEfficiency.C")

    ROOT.set_massRange(massLo, massHi, massBin)

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

    param_mass = ROOT.RooRealVar("m", "mass", massLo, massHi)


    mframe = ROOT.RooPlot("mframe","m frame", param_mass, massLo, massHi, massBin)

    # legend = ROOT.TLegend(1-margin_right-0.4, 0.5, 1-margin_right-0.01, 1-margin_top-0.01
    legend = ROOT.TLegend(margin_left+0.025, margin_bottom+0.03, 1-margin_right-0.075, 0.43)

    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    
    if True:
        model = "DasPlusExp"
        model = "Das"
        functions = []
        entries=[]
        for i, params in enumerate([
        # (90., 20, 1.5, 1.5, -0.1, 0.95),
        # (110., 20, 1.5, 1.5, -0.1, 0.95),
        # (90., 12, 1.5, 1.5, -0.1, 0.95),
        # (90., 20, 0.5, 0.5, -0.1, 0.95),
        # (90., 20, 1.5, 1.5, -0.1, 0.5),
        # (90., 20, 1.5, 1.5, -0.05, 0.5),
        # (130., 10, 1.0, 0.2),#, -0.3, 0.90),
        # (88., 17, 0.02, 0.6),#, -0.01, 0.95),
        # (130., 10, 2.0, 0.1),#, -0.1, 0.6),
        # (100., 10, 2.5, 0.1),#, -0.1, 0.85),
        # (100., 10, 8.0, 0.2),#, -0.01, 0.3),
        # (86., 36, 0.05, 3.0),#, -0.01, 0.6),
        (95., 30, 2.0, 0.5),#, -0.3, 0.90),
        (85., 25, 1.0, 0.6),#, -0.01, 0.95),
        (94., 50, 0.02, 0.5),#, -0.1, 0.6),
        (70., 20, 8., 0.2),#, -0.1, 0.85),
        (56., 40, 0.1, 1.2),#, -0.01, 0.3),
        (80., 30, 1.5, 1.0),#, -0.01, 0.6),
        ]):

            entry = ROOT.TLegendEntry()
            entry.SetTextSize(textsize)
            entry.SetLabel("#mu={0}, #sigma={1}".format(int(params[0]),params[1])
                +", k_{l}="+str(params[2])
                +", k_{h}="+str(params[3]))
            # +", t="+str(params[4])
            # +", f={0}".format(params[5]))
            entry.SetOption("L")
            entry.SetLineStyle(1)
            entry.SetLineColor(i+1)
            entry.SetLineWidth(2)
            entry.SetTextAlign(12)
            legend.GetListOfPrimitives().Add(entry)
            entries.append(entry)

            fX = ROOT.CDas(param_mass, 0, 0)
            fX.mean.setVal(params[0])
            fX.sigma.setVal(params[1])
            fX.kLo.setVal(params[2])
            fX.kHi.setVal(params[3])
            # fX.t1.setVal(params[4])
            # fX.frac.setVal(params[5])
            fX.model.plotOn(mframe)
            mframe.getAttLine().SetLineColor(i+1)
            functions.append(fX)

    if False:

        model = "CMSShape"
        functions = []
        entries=[]
        for i, params in enumerate([
            # (90., 0.05, 0.1),
            # (110., 0.05, 0.1),
            # (90., 0.05, 0.05),
            # (90., 0.05, 0.15),
            # (90., 0.0, 0.1),
            # (90., 0.1, 0.1),
            (66., 0.04, 0.015),
            (66., 0.04, 0.03),
            (46., 0.03, 0.025),
            (56., 0.01, 0.01),
            (56., 0.025, 0.03),
            (56., 0.05, 0.02),
        ]):
    
            entry = ROOT.TLegendEntry()
            entry.SetTextSize(textsize)
            entry.SetLabel("#alpha={0}, #beta={1}, #gamma={2}".format(int(params[0]), params[1], params[2]))
            entry.SetOption("L")
            entry.SetLineStyle(1)
            entry.SetLineColor(i+1)
            entry.SetLineWidth(2)
            entry.SetTextAlign(12)
            legend.GetListOfPrimitives().Add(entry)
            entries.append(entry)
            fX = ROOT.CRooCMSShape(param_mass, 0, 0, massLo, massHi)
            fX.alpha.setVal(params[0])
            fX.beta.setVal(params[1])
            fX.gamma.setVal(params[2])
            # fX.model.SetLineColor(i+1)
            fX.model.plotOn(mframe)
            mframe.getAttLine().SetLineColor(i+1)
            functions.append(fX)
    
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

    canvas.SaveAs(output+"/"+model+".png")
    canvas.SaveAs(output+"/"+model+".eps")
    canvas.Close()
