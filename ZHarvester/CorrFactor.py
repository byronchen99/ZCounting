###
#   Script to measure the correlation factor in data from the correlation between the HLT muons
#
###

import ROOT
import os
import numpy as np
import json
import pdb
import uncertainties as unc
from roofit.utils import load_makros

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))

# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

### Helper functions ###
# ---------------------->>>

# <<<------------------------


################################################################################
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input directory to the histogram file",
                        required=True, nargs="+")
    parser.add_argument("-o", "--output", help="Output directory to the histogram file",
                        default="./Test")
    parser.add_argument("-c", "--collect", help="Collect fits results only",
                        default=False, action="store_true",)
    args = parser.parse_args()

    input = args.input
    output = args.output

    signal_templates_path = "/nfs/dust/cms/user/dwalter/ZCounting/CMSSW_10_6_26/src/ZCounting/ZHarvester/res/sigTemplates/V13_02_02/"

    bkgModel = 6

    load_makros(bkgModel, 30, 2.4, 2017, 70, 250, 180, 0, 75)

    ROOT.set_massRange(70, 250, 180)

    ROOT.set_ptCut(30)

    if not os.path.isdir(output):
        os.mkdir(output)
    ROOT.set_luminosity(0)

    info = {}
    for eraMC, eraData, bins in [
        ("Summer16preVFP","2016preVFP", [0,13,18,75]),
        ("Summer16postVFP","2016postVFP", [0,15,20,75]),
        # ("Summer16preVFP","2016preVFP", [0,16,75]),
        # ("Summer16postVFP","2016postVFP", [0,75]),
        # ("Fall17", "2017", [0,18,23,28,35,75]),
        ("Fall17", "2017", [0,17,21,25,29,35,75]),
        # ("Fall17", "2017", [0,16,19,22,25,28,32,38,75]),
        ("Autumn18", "2018", [0,17,21,25,29,35,75]),
        # ("Autumn18", "2018", [0,22,30,75]),
    ]:
        info[eraData] = {}

        signal_template= signal_templates_path + "/ZSignalTemplate-V13_02_02-{0}-DYJetsToLL_M_50_NLO.root".format(eraMC)

        this_input = filter(lambda x: eraData in x, input)[0]

        file_ = ROOT.TFile(this_input,"READ")

        subdir = output +"/"+ eraData
        if not os.path.isdir(subdir):
            os.mkdir(subdir)
        ROOT.set_output(subdir)

        for region in (
            # "BB",
            # "BE",
            # "EE",
            "I",
            ):
            extension = "_HLT28"
            # extension = ""
            info[eraData][region] = {}

            for ibin in range(len(bins)-1):

                h0_ = file_.Get("hist_events{0}_0HLT_{1}_nPV_mass".format(extension, region)
                    ).ProjectionX("h1D_0HLT_{0}_mass".format(region),bins[ibin]+1,bins[ibin+1])
                h1_ = file_.Get("hist_events{0}_1HLT_{1}_nPV_mass".format(extension, region)
                    ).ProjectionX("h1D_1HLT_{0}_mass".format(region),bins[ibin]+1,bins[ibin+1])
                h2_ = file_.Get("hist_events{0}_2HLT_{1}_nPV_mass".format(extension, region)
                    ).ProjectionX("h1D_2HLT_{0}_mass".format(region),bins[ibin]+1,bins[ibin+1])

                h0_.SetDirectory(0)
                h1_.SetDirectory(0)
                h2_.SetDirectory(0)

                # load nPV distribution for center of mass
                hx_ = file_.Get("hist_events{0}_2HLT_{1}_nPV_mass".format(extension, region)
                    ).ProjectionY("h1D_2HLT_{0}_nPV".format(region))
                
                # file_.Get("hist_events_reco_{0}_nPV".format(region)).Clone("hist_events_reco_{0}_nPV_clone".format(region))
                for i in range(0,bins[ibin]):
                    hx_.SetBinContent(i, 0)
                    hx_.SetBinError(i, 0)
                for i in range(bins[ibin+1]+1,bins[-1]+1):
                    hx_.SetBinContent(i, 0)
                    hx_.SetBinError(i, 0)

                xCenter_ = hx_.GetMean()
                # xy_ = np.array([(i, hx_.GetBinContent(i)) for i in range(bins[ibin]+1,bins[ibin+1]+1)])
                # x_ = np.array([x[0] for x in xy_])
                # y_ = np.array([x[1] for x in xy_])
                #
                # xx_ = np.arange(x_[0],x_[-1],0.1)
                # yy_ = np.interp(xx_, x_, y_)
                # xLo_, xCenter_, xHi_ = np.percentile(yy_, [15.865, 50.0, 84.135], interpolation='nearest')
                # xLo_ = xx_[np.where(yy_==xLo_)[0]][0]
                # xCenter_ = xx_[np.where(yy_==xCenter_)[0]][0]
                # xHi_ = xx_[np.where(yy_==xHi_)[-1]][0]

                xLo_ = abs(bins[ibin]+1 - xCenter_)
                xHi_ = abs(bins[ibin+1] - xCenter_)

                if not args.collect:
                    ROOT.calculateHLTCorrelation(h0_, h1_, h2_, ibin, region, 2, bkgModel, 0, signal_template)

                # collect results and write in dictionary
                file = ROOT.TFile(subdir+"/workspace_{0}_{1}.root".format(region,ibin),"READ")
                workspace = file.Get("workspace")
                corr = unc.ufloat(workspace.var("c").getVal(), workspace.var("c").getError())

                bin = "{0}to{1}".format(bins[ibin]+1, bins[ibin+1])
                info[eraData][region][bin] = {
                    "c": workspace.var("c").getVal(),
                    "c_err": workspace.var("c").getError(),
                    "yLo": workspace.var("c").getErrorLo() ,
                    "yHi": workspace.var("c").getErrorHi(),
                    "x": xCenter_,
                    "xLo": xLo_,
                    "xHi": xHi_
                }

        file_.Close()

    # save result.json file
    with open(output+"/result.json","w") as outfile:
        json.dump(info, outfile, indent=4, sort_keys=True)
