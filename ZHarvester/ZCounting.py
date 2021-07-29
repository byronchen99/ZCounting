import logging as log
import ROOT
import pandas as pd
import glob
import os
import numpy as np
import json
import pdb
from uncertainties import ufloat
import threading, time

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime, getMCCorrection, unorm, pquad

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None
# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

### Helper functions ###
# ---------------------->>>

def get_correlationfactors(data):

    with open("res/correlationfactors.json") as file:
        res = json.load(file)

    retvals = []
    for region in ("BB","BE","EE"):
        factor = pquad(data['avgpu'], *res[era][region]["nominal"])
        factor = sum(factor * data['recorded(/pb)'].values) / sum(data['recorded(/pb)'].values)

        retvals.append(factor)

    return retvals


# <<<------------------------
def measurement(m,
    h2HLTBB_, h2HLTBE_, h2HLTEE_,
    h1HLTBB_, h1HLTBE_, h1HLTEE_,
    cBB_, cBE_, cEE_,
    h1SITBPass, h1SITBFail, h1SITEPass, h1SITEFail,
    h1TrkBPass, h1TrkBFail, h1TrkEPass, h1TrkEFail,
    h1StaBPass, h1StaBFail, h1StaEPass, h1StaEFail,
    sigmod_yield=1, bkgmod_yield=6,
    sigtempl_yield="", bkgtempl_yield="",
    sigmod_sel=[1,1,1,1], bkgmod_sel=[6,6,6,6],
    sigtempl_sel="", bkgshape_sel="",
    sigmod_trk=[1,1,1,1], bkgmod_trk=[6,6,6,6],
    sigtempl_trk="", bkgshape_trk="",
    sigmod_sta=[1,1,1,1], bkgmod_sta=[6,6,6,6],
    sigtempl_sta="", bkgshape_sta="",
):

    TrkeffresB_m = ROOT.calculateDataEfficiency(h1TrkBPass, h1TrkBFail,
                                                m, "Trk", 0, sigmod_trk[0], bkgmod_trk[0], sigmod_trk[1], bkgmod_trk[1],
                                                0, sigtempl_trk, bkgshape_trk)
    TrkeffresE_m = ROOT.calculateDataEfficiency(h1TrkEPass, h1TrkEFail,
                                                m, "Trk", 1, sigmod_trk[2], bkgmod_trk[2], sigmod_trk[3], bkgmod_trk[3],
                                                0, sigtempl_trk, bkgshape_trk)

    StaeffresB_m = ROOT.calculateDataEfficiency(h1StaBPass, h1StaBFail,
                                                m, "Sta", 0, sigmod_sta[0], bkgmod_sta[0], sigmod_sta[1], bkgmod_sta[1],
                                                0, sigtempl_sta, bkgshape_sta)
    StaeffresE_m = ROOT.calculateDataEfficiency(h1StaEPass, h1StaEFail,
                                                m, "Sta", 1, sigmod_sta[2], bkgmod_sta[2], sigmod_sta[3], bkgmod_sta[3],
                                                0, sigtempl_sta, bkgshape_sta)

    SeleffresB_m = ROOT.calculateDataEfficiency(h1SITBPass, h1SITBFail,
                                                m, "Sel", 0, sigmod_sel[0], bkgmod_sel[0], sigmod_sel[1], bkgmod_sel[1],
                                                0, sigtempl_sel, bkgshape_sel)
    SeleffresE_m = ROOT.calculateDataEfficiency(h1SITEPass, h1SITEFail,
                                                m, "Sel", 1, sigmod_sel[2], bkgmod_sel[2], sigmod_sel[3], bkgmod_sel[3],
                                                0, sigtempl_sel, bkgshape_sel)

    TrkeffB_m = ufloat(TrkeffresB_m[0], (TrkeffresB_m[1]+TrkeffresB_m[2])/2.)
    TrkeffE_m = ufloat(TrkeffresE_m[0], (TrkeffresE_m[1]+TrkeffresE_m[2])/2.)
    StaeffB_m = ufloat(StaeffresB_m[0], (StaeffresB_m[1]+StaeffresB_m[2])/2.)
    StaeffE_m = ufloat(StaeffresE_m[0], (StaeffresE_m[1]+StaeffresE_m[2])/2.)
    SeleffB_m = ufloat(SeleffresB_m[0], (SeleffresB_m[1]+SeleffresB_m[2])/2.)
    SeleffE_m = ufloat(SeleffresE_m[0], (SeleffresE_m[1]+SeleffresE_m[2])/2.)

    # ZtoMuMu efficiency purely from data
    ZBBeff = StaeffB_m * StaeffB_m * TrkeffB_m * TrkeffB_m * SeleffB_m * SeleffB_m
    ZBEeff = StaeffB_m * StaeffE_m * TrkeffB_m * TrkeffE_m * SeleffB_m * SeleffE_m
    ZEEeff = StaeffE_m * StaeffE_m * TrkeffE_m * TrkeffE_m * SeleffE_m * SeleffE_m

    res = {
        "SeleffB": SeleffB_m,
        "SeleffE": SeleffE_m,
        "TrkeffB": TrkeffB_m,
        "TrkeffE": TrkeffE_m,
        "StaeffB": StaeffB_m,
        "StaeffE": StaeffE_m,
        "SeleffB_chi2pass": SeleffresB_m[3],
        "SeleffB_chi2fail": SeleffresB_m[4],
        "SeleffE_chi2pass": SeleffresE_m[3],
        "SeleffE_chi2fail": SeleffresE_m[4],
        "TrkeffB_chi2pass": TrkeffresB_m[3],
        "TrkeffB_chi2fail": TrkeffresB_m[4],
        "TrkeffE_chi2pass": TrkeffresE_m[3],
        "TrkeffE_chi2fail": TrkeffresE_m[4],
        "StaeffB_chi2pass": StaeffresB_m[3],
        "StaeffB_chi2fail": StaeffresB_m[4],
        "StaeffE_chi2pass": StaeffresE_m[3],
        "StaeffE_chi2fail": StaeffresE_m[4],
        "ZBBeff": ZBBeff,
        "ZBEeff": ZBEeff,
        "ZEEeff": ZEEeff,
    }

    if True:
        # combined fit of BB, BE and EE in 1HLT and 2HLT category
        yields_m = ROOT.calculateHLTEfficiencyAndYield(
            h2HLTBB_, h1HLTBB_,
            h2HLTEE_, h1HLTEE_,
            h2HLTBE_, h1HLTBE_,
            m, sigmod_yield, bkgmod_yield, sigmod_yield, bkgmod_yield,
            cBB_, cBE_, cEE_,
            0, sigtempl_yield, bkgtempl_yield)

        res["HLTeffB"] = ufloat(yields_m[0], (yields_m[1]+yields_m[2])/2.)
        res["HLTeffE"] = ufloat(yields_m[3], (yields_m[4]+yields_m[5])/2.)

        res["zYieldBB"] = ufloat(yields_m[6], (yields_m[7] + yields_m[8])/2.)
        res["zYieldBE"] = ufloat(yields_m[9], (yields_m[10] + yields_m[11])/2.)
        res["zYieldEE"] = ufloat(yields_m[12], (yields_m[13] + yields_m[14])/2.)

        res["zYieldBB_chi2pass"] = yields_m[15]
        res["zYieldBB_chi2fail"] = yields_m[16]
        res["zYieldBE_chi2pass"] = yields_m[17]
        res["zYieldBE_chi2fail"] = yields_m[18]
        res["zYieldEE_chi2pass"] = yields_m[19]
        res["zYieldEE_chi2fail"] = yields_m[20]

    else:
        # Fit each region independently
        Hlt1BB_m = ROOT.getZyield(h1HLTBB_, m, sigmod_yield, bkgmod_yield, "BB", 0, sigtempl_yield)
        Hlt1BE_m = ROOT.getZyield(h1HLTBE_, m, sigmod_yield, bkgmod_yield, "BE", 0, sigtempl_yield)
        Hlt1EE_m = ROOT.getZyield(h1HLTEE_, m, sigmod_yield, bkgmod_yield, "EE", 0, sigtempl_yield)
        Hlt2BB_m = ROOT.getZyield(h2HLTBB_, m, sigmod_yield, bkgmod_yield, "BB", 1, sigtempl_yield)
        Hlt2BE_m = ROOT.getZyield(h2HLTBE_, m, sigmod_yield, bkgmod_yield, "BE", 1, sigtempl_yield)
        Hlt2EE_m = ROOT.getZyield(h2HLTEE_, m, sigmod_yield, bkgmod_yield, "EE", 1, sigtempl_yield)

        n1BB = ufloat(Hlt1BB_m[0], (Hlt1BB_m[1] + Hlt1BB_m[2])/2.)
        n1BE = ufloat(Hlt1BE_m[0], (Hlt1BE_m[1] + Hlt1BE_m[2])/2.)
        n1EE = ufloat(Hlt1EE_m[0], (Hlt1EE_m[1] + Hlt1EE_m[2])/2.)
        n2BB = ufloat(Hlt2BB_m[0], (Hlt2BB_m[1] + Hlt2BB_m[2])/2.)
        n2BE = ufloat(Hlt2BE_m[0], (Hlt2BE_m[1] + Hlt2BE_m[2])/2.)
        n2EE = ufloat(Hlt2EE_m[0], (Hlt2EE_m[1] + Hlt2EE_m[2])/2.)

        effHLT = lambda n1, n2, c: 1./c*(1./(1.+n1/n2))
        yieldZ = lambda n1, n2, c: n2*c*(1.+n1/(2.*n2))**2

        res["n1BB"] = n1BB
        res["n2BB"] = n2BB
        res["n1BE"] = n1BE
        res["n2BE"] = n2BE
        res["n1EE"] = n1EE
        res["n2EE"] = n2EE
        res["cBB"] = cBB_
        res["cBE"] = cBE_
        res["cEE"] = cEE_

        res["HLTeffB"] = effHLT(n1BB, n2BB, cBB_)
        res["HLTeffE"] = effHLT(n1EE, n2EE, cEE_)

        res["zYieldBB"] = yieldZ(n1BB, n2BB, cBB_)
        res["zYieldBE"] = yieldZ(n1BE, n2BE, cBE_)
        res["zYieldEE"] = yieldZ(n1EE, n2EE, cEE_)

        res["zYieldBB_chi2pass"] = Hlt2BB_m[3]
        res["zYieldBB_chi2fail"] = Hlt1BB_m[3]
        res["zYieldBE_chi2pass"] = Hlt2BE_m[3]
        res["zYieldBE_chi2fail"] = Hlt1BE_m[3]
        res["zYieldEE_chi2pass"] = Hlt2EE_m[3]
        res["zYieldEE_chi2fail"] = Hlt1EE_m[3]

    res["zDelBB"] = res["zYieldBB"]/res["ZBBeff"]
    res["zDelBE"] = res["zYieldBE"]/res["ZBEeff"]
    res["zDelEE"] = res["zYieldEE"]/res["ZEEeff"]

    return res

def ls_corrections(data_m, result_m, m, h2ZyieldBB, h2ZyieldBE, h2ZyieldEE):
    # --- per ls dataframe, one per measurement
    # for eta, hist in (("BB",h2ZyieldBB), ("BE",h2ZyieldBE), ("EE",h2ZyieldEE)):
        # if corr:
        #     # apply corrections
        #     data_m["eff{0}_mc".format(eta)] = result_m["Z{0}eff".format(eta)]#FIXME / corr['eff'+eta](data_m['avgpu'])
        # else:
        #     data_m["eff{0}_mc".format(eta)] = result_m["Z{0}eff".format(eta)]
        #
        # data_m['yield{0}'.format(eta)] = data_m['ls'].apply(lambda ls: hist.ProjectionY("", ls, ls, "e").Integral())
        #
        # data_m['zYield{0}'.format(eta)]  = data_m['yield{0}'.format(eta)] * result_m["zYield{0}_purity".format(eta)]
        # data_m['zDel{0}'.format(eta)]    = data_m['zYield{0}'.format(eta)] / result_m["Z{0}eff".format(eta)]
        # data_m['zDel{0}_mc'.format(eta)] = data_m['zYield{0}'.format(eta)] / data_m['eff{0}_mc'.format(eta)]


    data_m['time'] = data_m['time'].apply(lambda x: to_RootTime(x,currentYear))
    data_m['run'] = np.ones(len(data_m),dtype='int') * run
    data_m['fill'] = np.ones(len(data_m),dtype='int') * fill

    with open(outCSVDir + 'csvfile{0}_{1}.csv'.format(run, m), 'w') as file:
        data_m.to_csv(file, index=False)

    # --- per measurement result
    delLumi_m = data_m['delivered(/pb)'].sum()
    recLumi_m = data_m['recorded(/pb)'].sum()
    deadtime_m = recLumi_m / delLumi_m
    timeWindow_m = len(data_m) * secPerLS

    res = {
        "fill": fill,
        "run": run,
        "measurement": m,
        "tdate_begin": min(data_m['time']),
        "tdate_end": max(data_m['time']),
        # "yieldBB": data_m['yieldBB'].sum(),
        # "yieldBE": data_m['yieldBE'].sum(),
        # "yieldEE": data_m['yieldEE'].sum(),
        # "zYieldBB": data_m['zYieldBB'].sum(),
        # "zYieldBE": data_m['zYieldBE'].sum(),
        # "zYieldEE": data_m['zYieldEE'].sum(),
        # "zDelBB": data_m['zDelBB'].sum(),
        # "zDelBE": data_m['zDelBE'].sum(),
        # "zDelEE": data_m['zDelEE'].sum(),
        # "zDelBB_mc": data_m['zDelBB_mc'].sum(),
        # "zDelBE_mc": data_m['zDelBE_mc'].sum(),
        # "zDelEE_mc": data_m['zDelEE_mc'].sum(),
        # "xsecBB": data_m['zDelBB'].sum() / recLumi_m,
        # "xsecBE": data_m['zDelBE'].sum() / recLumi_m,
        # "xsecEE": data_m['zDelEE'].sum() / recLumi_m,
        # "xsecBB_mc": data_m['zDelBB_mc'].sum() / recLumi_m,
        # "xsecBE_mc": data_m['zDelBE_mc'].sum() / recLumi_m,
        # "xsecEE_mc": data_m['zDelEE_mc'].sum() / recLumi_m,
        "lumiDel": delLumi_m,
        "lumiRec": recLumi_m,
        "timewindow": timeWindow_m,
        "deadtime": deadtime_m,
        "pileUp": data_m['avgpu'].mean(),
        # "ZBBeff_mc": data_m['effBB_mc'].values.mean(),
        # "ZBEeff_mc": data_m['effBE_mc'].values.mean(),
        # "ZEEeff_mc": data_m['effEE_mc'].values.mean(),
    }
    return res

def load_histo(name_, fileName_, lumisections_=[0,], prefix_="", run_=0, suffix_=""):
    file_ = ROOT.TFile(fileName_)

    h_X_ls = file_.Get("{0}{1}".format(prefix_,name_)).Clone("{0}{1}{2}".format(prefix_,name_, suffix_))
    h_X_ls.SetDirectory(0)
    h_X = h_X_ls.ProjectionY("h_tmp_{0}_{1}_{2}".format(name_, run_, suffix_), lumisections_[0], lumisections_[0], "e")
    for ls in lumisections_[1:]:
        h_X.Add(h_X_ls.ProjectionY("h_tmp_{0}_{1}_{2}_{3}".format(name_, run_, ls, suffix_), ls, ls, "e"))

    # create new histogram in correct bin range
    hNew = ROOT.TH1D("h_mass_{0}_{1}_{2}_{3}".format(
        name_, run_, lumisections_[0], suffix_),
        "",MassBin_, MassMin_, MassMax_)

    for ibin in range(h_X.GetNbinsX()):
        binCenter = h_X.GetBinCenter(ibin)
        if binCenter < MassMin_:
            continue
        elif binCenter > MassMax_:
            break
        else:
            content = h_X.GetBinContent(ibin)
            newBin = hNew.FindBin(binCenter)
            hNew.SetBinContent(newBin, content)

    hNew.SetDirectory(0)

    return hNew


## Write big per measurement csv file
def writeSummaryCSV(outCSVDir):

    log.info("===Writing overall CSV file")
    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile??????.csv'))
    df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)

    with open(outCSVDir + '/Mergedcsvfile_perMeasurement.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

    # log.info("===Writing overall CSV file per LS")
    # rateFileList = sorted(glob.glob(outCSVDir + '/csvfile*_*.csv'))
    # # csvList = []
    # # for m in rateFileList:
    # #     csv = pd.read_csv(m)
    # #     measurement = int(m.split("_")[-1][:-4])
    # #     fill = int(m.split("_")[-2].split("csvfile")[-1])
    # #     csv["measurement"] = measurement
    # #     csvList.append(csv)
    # # df_merged = pd.concat(csvList, ignore_index=True)
    #
    # df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)
    #
    # with open(outCSVDir + '/Mergedcsvfile_perLS.csv', 'w') as file:
    #     df_merged.to_csv(file, index=False)

## Get root file with histograms
def getFileName(directory, run):
    # check if run was processed already
    log.info("===Loading input DQMIO.root file...")
    eosFileList = glob.glob(directory + '/*/*' + str(run) + '*root')
    if not len(eosFileList) > 0:
        log.info("The file does not yet exist for run: " + str(run))
        return None
    elif len(eosFileList) > 1:
        log.info("Multiple files found for run: " + str(run))
        return None
    else:
        return eosFileList[0]

################################################################################
if __name__ == '__main__':
    import argparse
    import os

    cmsswbase = os.environ['CMSSW_BASE']

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beginRun", help="first run to analyze [%(default)s]", type=int, default=272007)
    parser.add_argument("-e", "--endRun", help="analyze stops when comes to this run [%(default)s]", type=int, default=1000000)
    parser.add_argument('--mcCorrections', default="default", type=str,
                        help='specify .json file with MC corrections for muon correlations')
    parser.add_argument("-v", "--verbose", help="increase logging level from INFO to DEBUG", default=False,
                        action="store_true")
    parser.add_argument("-c", "--writeSummaryCSV", help="produce merged CSV with all runs", default=True)
    parser.add_argument("-i", "--dirDQM", help="Directory to the input root files from the DQM Offline module",
                        default="default")
    parser.add_argument("--byLsCSV", help="ByLs csv input generated by testBril.sh",
                        default="default")
    parser.add_argument("--sigTemplates", help="Directory to root file for Z mass template ('None': take analytic function) ",
                        default="default", type=str)
    parser.add_argument("--bkgTemplates", help="Directory to root file for bkg template ('None': take analytic function) ",
                        default="None", type=str)
    parser.add_argument('--ptCut', type=float, default=30.,
                        help='specify lower pt cut on tag and probe muons')
    parser.add_argument('--mass', nargs=3, metavar=('LOW', 'HIGH', 'NUMBER'), default=(70,250,180), type=int,
                        help='specify mass range for tag and probe muon pairs')
    parser.add_argument('--LumiPerMeasurement', default=20, type=float,
                        help='specify amount of luminosity per measurement in pb-1')
    parser.add_argument('--inclusive', default=False, action="store_true",
                        help='specify whether or not to do an inclusive fit of the specified runs')
    parser.add_argument("-o", "--dirOut", help="where to store the output files", default="./")

    args = parser.parse_args()
    if args.verbose:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
    else:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.INFO)


    ########################################
    # link to resouces
    prefix_dqm="" #"DQMData/Run {0}/ZCounting/Run summary/Histograms/".format(run)
    resPath = cmsswbase + "/src/ZCounting/ZHarvester/res/"
    if args.beginRun >= 272007 and args.beginRun < 278769:      # 2016 pre VFP
        byLsCSV         = resPath+"/FillByLs_RunII.csv"
        mcFile          = resPath+"/mcCorrections/corrections_nPU_2016preVFP.p"
        sigTemplates    = resPath+"/sigTemplates_V14_02/DYJetsToLL_M-50_2016preVFP"
        bkgTemplates    = resPath+"/qcdTemplates/Run2016preVFP.root"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11_02/Histograms_UL_RunII"
        era = "2016preVFP"
    elif args.beginRun < 294645:    # 2016 post VFP
        byLsCSV         = resPath+"/FillByLs_RunII.csv"
        mcFile          = resPath+"/mcCorrections/corrections_nPU_2016postVFP.p"
        sigTemplates    = resPath+"/sigTemplates_V14_02/DYJetsToLL_M-50_2016postVFP"
        bkgTemplates    = resPath+"/qcdTemplates/Run2016postVFP.root"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11_02/Histograms_UL_RunII"
        era = "2016postVFP"
    elif args.beginRun > 297020 and args.beginRun < 306828:     # 2017
        byLsCSV         = resPath+"/FillByLs_RunII.csv"
        mcFile          = resPath+"/mcCorrections/corrections_nPU_2017.p"
        sigTemplates    = resPath+"/sigTemplates_V14_02/DYJetsToLL_M-50_2017"
        bkgTemplates    = resPath+"/qcdTemplates/Run2017.root"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11_02/Histograms_UL_RunII"
        era = "2017"
    elif args.beginRun >= 306926 and args.beginRun < 307083:    # 2017 H
        byLsCSV         = resPath+"/FillByLs_2017_lowPU.csv"
        mcFile          = resPath+"/mcCorrections_SingleMu25/corrections_nPU_2017.p"
        sigTemplates    = resPath+"/sigTemplates_V14_02/DYJetsToLL_M-50_2017H"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11_02/Histograms_UL2017H_Pt30"
        era = "2017H"
    elif args.beginRun >= 315252:                               # 2018
        byLsCSV         = resPath+"/FillByLs_RunII.csv"
        mcFile          = resPath+"/mcCorrections/corrections_nPU_2018.p"
        sigTemplates    = resPath+"/sigTemplates_V14_02/DYJetsToLL_M-50_2018"
        bkgTemplates    = resPath+"/qcdTemplates/Run2018.root"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11_02/Histograms_UL_RunII"
        era = "2018"
    else:
        byLsCSV         = None
        mcFile          = None
        sigTemplates    = None
        eosDir          = None

    eosDir          = eosDir        if args.dirDQM        == "default"   else args.dirDQM
    byLsCSV         = byLsCSV       if args.byLsCSV       == "default"   else args.byLsCSV
    sigTemplates    = sigTemplates  if args.sigTemplates  == "default"   else args.sigTemplates
    bkgTemplates    = bkgTemplates  if args.bkgTemplates  == "default"   else args.bkgTemplates
    mcFile          = mcFile        if args.mcCorrections == "default"   else args.mcCorrections

    print("----------------------------------")
    print("Use eosDir:       {0}".format(eosDir))
    print("Use byLsCSV:      {0}".format(byLsCSV))
    print("Use sigTemplates: {0}".format(sigTemplates))
    print("Use bkgTemplates: {0}".format(bkgTemplates))
    print("Use mcFile:       {0}".format(mcFile))
    print("Mass range from:  {0} to {1}".format(*args.mass))
    print("Lumi per Measurement:  {0}".format(args.LumiPerMeasurement))
    print("----------------------------------")

    ########## Input configuration ##########
    # ByLS csv inputs generated by testBRIL.sh
    byLS_filelist = glob.glob(byLsCSV)
    byLS_filelist.sort(key=os.path.getmtime)
    byLS_filename = byLS_filelist[-1]
    log.info("The brilcalc csv file: " + str(byLS_filename))

    # inputs: to build MC*Gaussian template for efficiency fitting
    fit_options={}
    if sigTemplates and sigTemplates != "None":
        _optns = {
            'sigmod_yield':2,
            'sigtempl_yield':sigTemplates+"/template_ZYield.root",
            # 'sigmod_hlt':[2,2,2,2],
            # 'sigtempl_hlt':sigTemplates+"/template_HLT.root",
            'sigmod_sel':[2,2,2,2],
            'sigtempl_sel':sigTemplates+"/template_Sel.root",
            'sigmod_trk':[2,2,2,2],
            'sigtempl_trk':sigTemplates+"/template_Trk.root",
            'sigmod_sta':[2,2,2,2],
            'sigtempl_sta':sigTemplates+"/template_Sta.root",
            }
        fit_options.update(_optns)

    if bkgTemplates and bkgTemplates != "None":

        _optns = {
            # "bkgmod_hlt": [7,7,7,7],
            # "bkgshape_hlt": bkgTemplates,
            "bkgmod_sel": [7,7,7,7],
            "bkgshape_sel": bkgTemplates,
            "bkgmod_trk": [7,7,7,7],
            "bkgshape_trk": bkgTemplates,
            "bkgmod_sta": [7,7,7,7],
            "bkgshape_sta": bkgTemplates,
            }
        fit_options.update(_optns)

    outDir = args.dirOut if args.dirOut.endswith("/") else args.dirOut+"/"
    if not os.path.isdir(outDir):
        os.mkdir(outDir)

    outCSVDir = outDir+"csvFiles/"
    if not os.path.isdir(outCSVDir):
        try:
            os.mkdir(outCSVDir)
        except OSError:
            log.warning("directory already exists ...")

    # if mcFile:
    #     print("Use MC corrections from: ")
    #     print("file: "+mcFile)
    #     corr = getMCCorrection(mcFile)
    # else:
    #     corr = None

    ########### Constant settings ##########
    secPerLS = float(23.3)
    currentYear = 2017

    maximumLS = 5000
    # LSperMeasurement = 100  # required number of lumi sections per measurement
    LumiPerMeasurement = args.LumiPerMeasurement  # minimum recorded lumi for one measurement in pb-1

    #nonfiguration of fit models

    log.info("Loading C marco...")
    ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(
        __file__)) + "/calculateDataEfficiency.C")  # load function getZyield(...) and calculateDataEfficiency(...)

    # if args.bkgSubtract:
    #     MassBin_ = 3
    #     MassMin_ = -1.5
    #     MassMax_ = 1.5
    # else:

    MassMin_ = int(args.mass[0])
    MassMax_ = int(args.mass[1])
    MassBin_ = int(args.mass[2])

    ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

    ROOT.set_ptCut(args.ptCut)

    log.info("Loading input byls csv...")
    byLS_file = open(str(byLS_filename))
    byLS_lines = byLS_file.readlines()
    byLS_data = pd.read_csv(byLS_filename, sep=',', low_memory=False,
                       skiprows=lambda x: byLS_lines[x].startswith('#') and not byLS_lines[x].startswith('#run'))
    log.debug("%s", byLS_data.axes)

    log.info("formatting csv file...")    # formatting the csv
    byLS_data[['run', 'fill']] = byLS_data['#run:fill'].str.split(':', expand=True).apply(pd.to_numeric)
    byLS_data['ls'] = byLS_data['ls'].str.split(':', expand=True)[0].apply(pd.to_numeric)
    byLS_data = byLS_data.drop(['#run:fill', 'hltpath', 'source'], axis=1)

    if 'delivered(/ub)' in byLS_data.columns.tolist():  # convert to /pb
        byLS_data['delivered(/ub)'] = byLS_data['delivered(/ub)'].apply(lambda x: x / 1000000.)
        byLS_data['recorded(/ub)'] = byLS_data['recorded(/ub)'].apply(lambda x: x / 1000000.)
        byLS_data = byLS_data.rename(index=str, columns={'delivered(/ub)': 'delivered(/pb)', 'recorded(/ub)': 'recorded(/pb)'})
    elif 'delivered(/fb)' in byLS_data.columns.tolist():  # convert to /pb
        byLS_data['delivered(/fb)'] = byLS_data['delivered(/fb)'].apply(lambda x: x * 1000.)
        byLS_data['recorded(/fb)'] = byLS_data['recorded(/fb)'].apply(lambda x: x * 1000.)
        byLS_data = byLS_data.rename(index=str, columns={'delivered(/fb)': 'delivered(/pb)', 'recorded(/fb)': 'recorded(/pb)'})

    # if there are multiple entries of the same ls (for example from different triggers), only keep the one with the highest luminosity.
    byLS_data = byLS_data.sort_values(['fill', 'run', 'ls', 'delivered(/pb)', 'recorded(/pb)'])
    byLS_data = byLS_data.drop_duplicates(['fill', 'run', 'ls'])

    #####################################

    if args.inclusive:
        #perform one inclusive fit for the full data
        recLumi = 0

        # For fitting
        h2HLTBB = ROOT.TH1D("h_mass_2HLTBB_Z","",MassBin_, MassMin_, MassMax_)
        h2HLTBE = ROOT.TH1D("h_mass_2HLTBE_Z","",MassBin_, MassMin_, MassMax_)
        h2HLTEE = ROOT.TH1D("h_mass_2HLTEE_Z","",MassBin_, MassMin_, MassMax_)
        h1HLTBB = ROOT.TH1D("h_mass_1HLTBB_Z","",MassBin_, MassMin_, MassMax_)
        h1HLTBE = ROOT.TH1D("h_mass_1HLTBE_Z","",MassBin_, MassMin_, MassMax_)
        h1HLTEE = ROOT.TH1D("h_mass_1HLTEE_Z","",MassBin_, MassMin_, MassMax_)

        hSITpassB = ROOT.TH1D("h_mass_SIT_pass_central","",MassBin_, MassMin_, MassMax_)
        hSITfailB = ROOT.TH1D("h_mass_SIT_fail_central","",MassBin_, MassMin_, MassMax_)
        hSITpassE = ROOT.TH1D("h_mass_SIT_pass_forward","",MassBin_, MassMin_, MassMax_)
        hSITfailE = ROOT.TH1D("h_mass_SIT_fail_forward","",MassBin_, MassMin_, MassMax_)
        hTrkpassB = ROOT.TH1D("h_mass_Trk_pass_central","",MassBin_, MassMin_, MassMax_)
        hTrkfailB = ROOT.TH1D("h_mass_Trk_fail_central","",MassBin_, MassMin_, MassMax_)
        hTrkpassE = ROOT.TH1D("h_mass_Trk_pass_forward","",MassBin_, MassMin_, MassMax_)
        hTrkfailE = ROOT.TH1D("h_mass_Trk_fail_forward","",MassBin_, MassMin_, MassMax_)
        hStapassB = ROOT.TH1D("h_mass_Sta_pass_central","",MassBin_, MassMin_, MassMax_)
        hStafailB = ROOT.TH1D("h_mass_Sta_fail_central","",MassBin_, MassMin_, MassMax_)
        hStapassE = ROOT.TH1D("h_mass_Sta_pass_forward","",MassBin_, MassMin_, MassMax_)
        hStafailE = ROOT.TH1D("h_mass_Sta_fail_forward","",MassBin_, MassMin_, MassMax_)

        for run, data_run in byLS_data.groupby('run'):
            if run < int(args.beginRun) or run >= int(args.endRun):
                continue

            LSlist = data_run.query('ls <= {0}'.format(maximumLS))['ls'].values.tolist()
            recLumi += sum(data_run.loc[data_run['ls'].isin(LSlist)]['recorded(/pb)'].values)

            eosFile = getFileName(eosDir, run)      # >>> histograms binned in mass
            if eosFile is None:
                log.info("Continue")
                continue

            def load(name_):
                return load_histo(name_, eosFile, LSlist, prefix_dqm, run, "new")

            h2HLTBB.Add(load("h_mass_2HLTBB_Z"))
            h2HLTBE.Add(load("h_mass_2HLTBE_Z"))
            h2HLTEE.Add(load("h_mass_2HLTEE_Z"))
            h1HLTBB.Add(load("h_mass_1HLTBB_Z"))
            h1HLTBE.Add(load("h_mass_1HLTBE_Z"))
            h1HLTEE.Add(load("h_mass_1HLTEE_Z"))

            hSITpassB.Add(load("h_mass_SIT_pass_central"))
            hSITfailB.Add(load("h_mass_SIT_fail_central"))
            hSITpassE.Add(load("h_mass_SIT_pass_forward"))
            hSITfailE.Add(load("h_mass_SIT_fail_forward"))
            hTrkpassB.Add(load("h_mass_Trk_pass_central"))
            hTrkfailB.Add(load("h_mass_Trk_fail_central"))
            hTrkpassE.Add(load("h_mass_Trk_pass_forward"))
            hTrkfailE.Add(load("h_mass_Trk_fail_forward"))
            hStapassB.Add(load("h_mass_Sta_pass_central"))
            hStafailB.Add(load("h_mass_Sta_fail_central"))
            hStapassE.Add(load("h_mass_Sta_pass_forward"))
            hStafailE.Add(load("h_mass_Sta_fail_forward"))


        outSubDir = outDir + "Fits_Inclusive/"
        os.mkdir(outSubDir)
        ROOT.set_output(outSubDir)
        ROOT.set_luminosity(recLumi)

        df = data_run.loc[data_run['ls'].isin(LSlist)]
        # get correlation factor for second muon
        cBB, cBE, cEE = get_correlationfactors(df)

        result = measurement(0,
            h2HLTBB, h2HLTBE, h2HLTEE,
            h1HLTBB, h1HLTBE, h1HLTEE,
            cBB, cBE, cEE,
            hSITpassB, hSITfailB, hSITpassE, hSITfailE,
            hTrkpassB, hTrkfailB, hTrkpassE, hTrkfailE,
            hStapassB, hStafailB, hStapassE, hStafailE,
            **fit_options
            )

    log.info("===Looping over runs... {0} to {1}".format(int(args.beginRun), int(args.endRun)))
    for run, data_run in byLS_data.groupby('run'):
        if run < int(args.beginRun) or run >= int(args.endRun):
            continue

        fill = data_run.drop_duplicates('fill')['fill'].values[0]
        LSlist = data_run.query('ls <= {0}'.format(maximumLS))['ls'].values.tolist()

        if not args.inclusive:
            # check if run was processed already
            outSubDir = outDir + "Run{0}/".format(run)
            if os.path.isdir(outSubDir):
                log.info("Run %i was already processed, skipping and going to next run", run)
                continue
            os.mkdir(outSubDir)

        ROOT.set_output(outSubDir)

        log.info("===Running Run %i", run)
        log.info("===Running Fill %i", fill)

        eosFile = getFileName(eosDir, run)      # >>> histograms binned in mass
        # eosFileSS = getFileName(eosDirSS, run)  # >>> histograms binned in charge of muon pair
        if eosFile is None:
            log.info("Continue")
            continue

        log.info("===Looping over measurements...")
        results = []
        while len(LSlist) > 0:  # begin next measurement "m"

            # merge data to one measuement if remaining luminosity is too less for two measuements
            mergeMeasurements = sum(data_run.loc[data_run['ls'].isin(LSlist)]['recorded(/pb)'].values) < 1.5 * LumiPerMeasurement

            # produce goodLSlist with ls that are used for one measurement
            goodLSlist = []
            recLumi = 0
            while len(LSlist) > 0:
                goodLSlist.append(LSlist[0])
                recLumi += (data_run[data_run['ls'] == LSlist[0]]['recorded(/pb)'].values)[0]
                del LSlist[0]
                 # if we have enough collected lumisections
                if not mergeMeasurements and recLumi >= LumiPerMeasurement:# or len(goodLSlist) >= LSperMeasurement):
                    break


            df = data_run.loc[data_run['ls'].isin(goodLSlist)]

            # get correlation factor for second muon
            cBB, cBE, cEE = get_correlationfactors(df)

            ### load histograms
            log.debug("Openning DQMIO.root file: %s", eosFile)

            def load(name_):
                return load_histo(name_, eosFile, goodLSlist, prefix_dqm, run, "new")

            if not args.inclusive:

                h2HLTBB = load("h_mass_2HLTBB_Z")
                h2HLTBE = load("h_mass_2HLTBE_Z")
                h2HLTEE = load("h_mass_2HLTEE_Z")
                h1HLTBB = load("h_mass_1HLTBB_Z")
                h1HLTBE = load("h_mass_1HLTBE_Z")
                h1HLTEE = load("h_mass_1HLTEE_Z")

                hSITpassB = load("h_mass_SIT_pass_central")
                hSITfailB = load("h_mass_SIT_fail_central")
                hSITpassE = load("h_mass_SIT_pass_forward")
                hSITfailE = load("h_mass_SIT_fail_forward")
                hTrkpassB = load("h_mass_Trk_pass_central")
                hTrkfailB = load("h_mass_Trk_fail_central")
                hTrkpassE = load("h_mass_Trk_pass_forward")
                hTrkfailE = load("h_mass_Trk_fail_forward")
                hStapassB = load("h_mass_Sta_pass_central")
                hStafailB = load("h_mass_Sta_fail_central")
                hStapassE = load("h_mass_Sta_pass_forward")
                hStafailE = load("h_mass_Sta_fail_forward")


                ROOT.set_output(outSubDir)
                ROOT.set_luminosity(recLumi)

                result = measurement(len(results),
                    h2HLTBB, h2HLTBE, h2HLTEE,
                    h1HLTBB, h1HLTBE, h1HLTEE,
                    cBB, cBE, cEE,
                    hSITpassB, hSITfailB, hSITpassE, hSITfailE,
                    hTrkpassB, hTrkfailB, hTrkpassE, hTrkfailE,
                    hStapassB, hStafailB, hStapassE, hStafailE,
                    **fit_options
                    )

            if result is None:
                continue

            result.update(ls_corrections(df, result, len(results), None, None, None))#hist_ZyieldBB, hist_ZyieldBE, hist_ZyieldEE))

            results.append(result)

        ## Write per measurement csv file - one per run
        log.info("===Writing per Run CSV file")
        results = pd.concat([pd.DataFrame([result]) for result in results])

        with open(outCSVDir + '/csvfile{0}.csv'.format(run), 'w') as file:
            results.to_csv(file, index=False)


    if args.writeSummaryCSV:
        writeSummaryCSV(outCSVDir)

    log.info("===Done")
