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
from ZUtils.python.utils import to_RootTime, getMCCorrection

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None
# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

### Helper functions ###
# ---------------------->>>

def unorm(x):
    # for counting experiments: define ufloat with poisson uncertainty
    return ufloat(x, np.sqrt(abs(x)))

# compute Z yield by counting
def count_bkg(histPos, histNeg):
    mm = unorm(histNeg.Integral())
    pp = unorm(histPos.Integral())

    # background are transfered via this formula to opposite sign region
    return 2 * (mm * pp)**0.5

def count_signal(histOS, histPos, histNeg):
    os = unorm(histOS.Integral())

    # transfer same sign events to opposite sign
    return os - count_bkg(histPos, histNeg)

# old implementation
# def count_signal(hist):
#     os = unorm(hist.Integral())
#
#     # transfer same sign events to opposite sign
#     return os - count_bkg(hist)
#
# def count_efficiency(hPass, hFail):
#     osPass = unorm(hPass.Integral(2,2))
#     ssPass = count_bkg(hPass)
#     osFail = unorm(hFail.Integral(2,2))
#     ssFail = count_bkg(hFail)


def count_purity(histOS, histPos, histNeg):
    os = unorm(histOS.Integral())   # opposite sign events -> signal + backrgound
    ss = count_bkg(histPos, histNeg)            # same sign events -> background
    return (os - ss) / os

    return (osPass - ssPass) / (osPass - ssPass + osFail - ssFail)
# <<<------------------------


def measurement(m,
    h1RecoBB, h1RecoBE, h1RecoEE,
    h1RecoBB_neg, h1RecoBE_neg, h1RecoEE_neg,
    h1RecoBB_pos, h1RecoBE_pos, h1RecoEE_pos,
    h1HLTBPass, h1HLTBFail, h1HLTEPass, h1HLTEFail,
    h1SITBPass, h1SITBFail, h1SITEPass, h1SITEFail,
    h1TrkBPass, h1TrkBFail, h1TrkEPass, h1TrkEFail,
    h1StaBPass, h1StaBFail, h1StaEPass, h1StaEFail,
    sigmod_yield=1, bkgmod_yield=6,
    sigtempl_yield="", bkghist_yield=0,
    sigmod_hlt=[1,1,1,1], bkgmod_hlt=[6,6,6,6],
    sigtempl_hlt="", bkgshape_hlt="",
    sigmod_sel=[1,1,1,1], bkgmod_sel=[6,6,6,6],
    sigtempl_sel="", bkgshape_sel="",
    sigmod_trk=[1,1,1,1], bkgmod_trk=[6,6,6,6],
    sigtempl_trk="", bkgshape_trk="",
    sigmod_sta=[1,1,1,1], bkgmod_sta=[6,6,6,6],
    sigtempl_sta="", bkgshape_sta="",
    ):

    # # compute Z yield from fits
    # # --->>>
    # ZyieldresBB_m = ROOT.getZyield(h1RecoBB, m, sigmod_yield, bkgmod_yield,
    #     "BB", sigtempl_yield, bkghist_yield)
    #
    # ZyieldresBE_m = ROOT.getZyield(h1RecoBE, m, sigmod_yield, bkgmod_yield,
    #     "BE", sigtempl_yield, bkghist_yield)
    #
    # ZyieldresEE_m = ROOT.getZyield(h1RecoEE, m, sigmod_yield, bkgmod_yield,
    #     "EE", sigtempl_yield, bkghist_yield)
    #
    # # check if there is was enough statistics
    # if np.isnan(ZyieldresBB_m[0]) or np.isnan(ZyieldresBE_m[0]) or np.isnan(ZyieldresEE_m[0]) or ZyieldresBB_m[0] * ZyieldresBE_m[0] * ZyieldresEE_m[0] == 0:
    #     return None
    #
    # ZyieldBB_m = ufloat(ZyieldresBB_m[0], (ZyieldresBB_m[1]+ZyieldresBB_m[2])/2.)
    # ZyieldBE_m = ufloat(ZyieldresBE_m[0], (ZyieldresBE_m[1]+ZyieldresBE_m[2])/2.)
    # ZyieldEE_m = ufloat(ZyieldresEE_m[0], (ZyieldresEE_m[1]+ZyieldresEE_m[2])/2.)
    #
    # ZpurityBB_m = ufloat(ZyieldresBB_m[4], (ZyieldresBB_m[5]+ZyieldresBB_m[6])/2.)
    # ZpurityBE_m = ufloat(ZyieldresBE_m[4], (ZyieldresBE_m[5]+ZyieldresBE_m[6])/2.)
    # ZpurityEE_m = ufloat(ZyieldresEE_m[4], (ZyieldresEE_m[5]+ZyieldresEE_m[6])/2.)
    # # <<<---

    # compute Z yield from counting events and subtract same sign events
    # --->>>
    ZyieldBB_m = count_signal(h1RecoBB, h1RecoBB_pos, h1RecoBB_neg)
    ZyieldBE_m = count_signal(h1RecoBE, h1RecoBE_pos, h1RecoBE_neg)
    ZyieldEE_m = count_signal(h1RecoEE, h1RecoEE_pos, h1RecoEE_neg)

    ZpurityBB_m = count_purity(h1RecoBB, h1RecoBB_pos, h1RecoBB_neg)
    ZpurityBE_m = count_purity(h1RecoBE, h1RecoBE_pos, h1RecoBE_neg)
    ZpurityEE_m = count_purity(h1RecoEE, h1RecoEE_pos, h1RecoEE_neg)
    # <<<---

    ### compute muon efficiencies
    HLTeffresB_m = ROOT.calculateDataEfficiency(h1HLTBPass, h1HLTBFail,
                                                m, "HLT", 0, sigmod_hlt[0], bkgmod_hlt[0], sigmod_hlt[1], bkgmod_hlt[1],
                                                0, sigtempl_hlt, bkgshape_hlt)
    HLTeffresE_m = ROOT.calculateDataEfficiency(h1HLTEPass, h1HLTEFail,
                                                m, "HLT", 1, sigmod_hlt[2], bkgmod_hlt[2], sigmod_hlt[3], bkgmod_hlt[3],
                                                0, sigtempl_hlt, bkgshape_hlt)

    SeleffresB_m = ROOT.calculateDataEfficiency(h1SITBPass, h1SITBFail,
                                                m, "Sel", 0, sigmod_sel[0], bkgmod_sel[0], sigmod_sel[1], bkgmod_sel[1],
                                                0, sigtempl_sel, bkgshape_sel)
    SeleffresE_m = ROOT.calculateDataEfficiency(h1SITEPass, h1SITEFail,
                                                m, "Sel", 1, sigmod_sel[2], bkgmod_sel[2], sigmod_sel[3], bkgmod_sel[3],
                                                0, sigtempl_sel, bkgshape_sel)

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


    HLTeffB_m = ufloat(HLTeffresB_m[0], (HLTeffresB_m[1]+HLTeffresB_m[2])/2.)
    HLTeffE_m = ufloat(HLTeffresE_m[0], (HLTeffresE_m[1]+HLTeffresE_m[2])/2.)
    SeleffB_m = ufloat(SeleffresB_m[0], (SeleffresB_m[1]+SeleffresB_m[2])/2.)
    SeleffE_m = ufloat(SeleffresE_m[0], (SeleffresE_m[1]+SeleffresE_m[2])/2.)
    TrkeffB_m = ufloat(TrkeffresB_m[0], (TrkeffresB_m[1]+TrkeffresB_m[2])/2.)
    TrkeffE_m = ufloat(TrkeffresE_m[0], (TrkeffresE_m[1]+TrkeffresE_m[2])/2.)
    StaeffB_m = ufloat(StaeffresB_m[0], (StaeffresB_m[1]+StaeffresB_m[2])/2.)
    StaeffE_m = ufloat(StaeffresE_m[0], (StaeffresE_m[1]+StaeffresE_m[2])/2.)

    # ZtoMuMu efficiency purely from data
    ZBBeff = (StaeffB_m * StaeffB_m * TrkeffB_m * TrkeffB_m * SeleffB_m * SeleffB_m * (1 - (1 - HLTeffB_m) * (1 - HLTeffB_m)))
    ZBEeff = (StaeffB_m * StaeffE_m * TrkeffB_m * TrkeffE_m * SeleffB_m * SeleffE_m * (1 - (1 - HLTeffB_m) * (1 - HLTeffE_m)))
    ZEEeff = (StaeffE_m * StaeffE_m * TrkeffE_m * TrkeffE_m * SeleffE_m * SeleffE_m * (1 - (1 - HLTeffE_m) * (1 - HLTeffE_m)))

    ZBB_m = ZyieldBB_m/ZBBeff
    ZBE_m = ZyieldBE_m/ZBEeff
    ZEE_m = ZyieldEE_m/ZEEeff

    res = {
        "zYieldBB": ZyieldBB_m,
        # "zYieldBB_chi2": ZyieldresBB_m[3],
        "zYieldBB_purity": ZpurityBB_m,
        "zYieldBE": ZyieldBE_m,
        # "zYieldBE_chi2": ZyieldresBE_m[3],
        "zYieldBE_purity": ZpurityBE_m,
        "zYieldEE": ZyieldEE_m,
        # "zYieldEE_chi2": ZyieldresEE_m[3],
        "zYieldEE_purity": ZpurityEE_m,
        "HLTeffB": HLTeffB_m,
        "HLTeffE": HLTeffE_m,
        "SeleffB": SeleffB_m,
        "SeleffE": SeleffE_m,
        "TrkeffB": TrkeffB_m,
        "TrkeffE": TrkeffE_m,
        "StaeffB": StaeffB_m,
        "StaeffE": StaeffE_m,
        "HLTeffB_chi2pass": HLTeffresB_m[3],
        "HLTeffB_chi2fail": HLTeffresB_m[4],
        "HLTeffE_chi2pass": HLTeffresE_m[3],
        "HLTeffE_chi2fail": HLTeffresE_m[4],
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
        # "ZBBeff_stat": ZBBeff.s,
        # "ZBEeff_stat": ZBEeff.s,
        # "ZEEeff_stat": ZEEeff.s,
        "zBB": ZBB_m,
        "zBE": ZBB_m,
        "zEE": ZBB_m,
        # "zBB_stat": ZBB_m.s,
        # "zBE_stat": ZBB_m.s,
        # "zEE_stat": ZBB_m.s,
    }

    return res

def ls_corrections(data_m, result_m, m, h2ZyieldBB, h2ZyieldBE, h2ZyieldEE):
    # --- per ls dataframe, one per measurement
    for eta, hist in (("BB",h2ZyieldBB), ("BE",h2ZyieldBE), ("EE",h2ZyieldEE)):

        if corr:
            data_m["eff{0}_mc".format(eta)] = result_m["Z{0}eff".format(eta)] / corr['eff'+eta](data_m['avgpu'])
        else:
            data_m["eff{0}_mc".format(eta)] = result_m["Z{0}eff".format(eta)]

        data_m['yield{0}'.format(eta)] = data_m['ls'].apply(lambda ls: hist.ProjectionY("", ls, ls, "e").Integral())

        data_m['zYield{0}'.format(eta)]  = data_m['yield{0}'.format(eta)] * result_m["zYield{0}_purity".format(eta)]
        data_m['zDel{0}'.format(eta)]    = data_m['zYield{0}'.format(eta)] / result_m["Z{0}eff".format(eta)]
        data_m['zDel{0}_mc'.format(eta)] = data_m['zYield{0}'.format(eta)] / data_m['eff{0}_mc'.format(eta)]

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
        "tdate_begin": min(data_m['time']),
        "tdate_end": max(data_m['time']),
        "yieldBB": data_m['yieldBB'].sum(),
        "yieldBE": data_m['yieldBE'].sum(),
        "yieldEE": data_m['yieldEE'].sum(),
        "zYieldBB": data_m['zYieldBB'].sum(),
        "zYieldBE": data_m['zYieldBE'].sum(),
        "zYieldEE": data_m['zYieldEE'].sum(),
        "zDelBB": data_m['zDelBB'].sum(),
        "zDelBE": data_m['zDelBE'].sum(),
        "zDelEE": data_m['zDelEE'].sum(),
        "zDelBB_mc": data_m['zDelBB_mc'].sum(),
        "zDelBE_mc": data_m['zDelBE_mc'].sum(),
        "zDelEE_mc": data_m['zDelEE_mc'].sum(),
        "xsecBB": data_m['zDelBB'].sum() / recLumi_m,
        "xsecBE": data_m['zDelBE'].sum() / recLumi_m,
        "xsecEE": data_m['zDelEE'].sum() / recLumi_m,
        "xsecBB_mc": data_m['zDelBB_mc'].sum() / recLumi_m,
        "xsecBE_mc": data_m['zDelBE_mc'].sum() / recLumi_m,
        "xsecEE_mc": data_m['zDelEE_mc'].sum() / recLumi_m,
        "lumiDel": delLumi_m,
        "lumiRec": recLumi_m,
        "timewindow": timeWindow_m,
        "deadtime": deadtime_m,
        "pileUp": data_m['avgpu'].mean(),
        "ZBBeff_mc": data_m['effBB_mc'].values.mean(),
        "ZBEeff_mc": data_m['effBE_mc'].values.mean(),
        "ZEEeff_mc": data_m['effEE_mc'].values.mean(),
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
    hNew = ROOT.TH1D("h_mass_{0}_{1}_{2}_{3}".format(name_, run_, ls, suffix_),"",MassBin_, MassMin_, MassMax_)

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

    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile*_*.csv'))
    df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)

    with open(outCSVDir + '/Mergedcsvfile_perLS.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

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
    parser.add_argument('--mcCorrections', default=None, type=str,
                        help='specify .json file with MC corrections for muon correlations')
    parser.add_argument("-v", "--verbose", help="increase logging level from INFO to DEBUG", default=False,
                        action="store_true")
    parser.add_argument("-c", "--writeSummaryCSV", help="produce merged CSV with all runs", default=True)
    parser.add_argument("-i", "--dirDQM", help="Directory to the input root files from the DQM Offline module",
                        default=None)
    parser.add_argument("--byLsCSV", help="ByLs csv input generated by testBril.sh",
                        default="/nfs/dust/cms/user/dwalter/data/Lumi/For_2017/FillByLs_2017.csv")
    parser.add_argument("--sigTemplates", help="Directory to root file for Z mass template",
                        default=None, type=str)
    parser.add_argument("--bkgTemplates", help="Directory to root file for bkg template",
                        default=None, type=str)
    parser.add_argument("--bkgSubtract", help="Directory to root files for bkg subtraction (Same sign (+) and (-) histograms)",
                        default=False, action="store_true")
    parser.add_argument('--ptCut', type=float, default=30.,
                        help='specify lower pt cut on tag and probe muons')
    parser.add_argument('--mass', nargs=2, metavar=('LOW', 'HIGH'), default=(56,116),
                        help='specify lower pt cut on tag and probe muons')
    parser.add_argument('--inclusive', default=False, action="store_true",
                        help='specify whether or not to do an inclusive fit of the specified runs')
    parser.add_argument("-o", "--dirOut", help="where to store the output files", default="./")

    args = parser.parse_args()
    if args.verbose:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
    else:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.INFO)

    prefix_dqm="" #"DQMData/Run {0}/ZCounting/Run summary/Histograms/".format(run)
    resPath = cmsswbase + "/src/ZCounting/ZHarvester/res/"

    ########################################
    if args.beginRun > 272007 and args.beginRun < 278769:      # 2016 pre VFP
        byLsCSV         = resPath+"/FillByLs_2016.csv"
        mcFile          = resPath+"/mcCorrections/corrections_nPU_2016preVFP.p"
        sigTemplates    = resPath+"/sigTemplates/DYJetsToLL_M-50_2016preVFP"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2016_Pt30M56to116"
        # eosDirSS        = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2016_Pt30M56to116_SS"
    elif args.beginRun < 294645:    # 2016 post VFP
        byLsCSV         = resPath+"/FillByLs_2016.csv"
        mcFile          = resPath+"/mcCorrections/corrections_nPU_2016postVFP.p"
        sigTemplates    = resPath+"/sigTemplates/DYJetsToLL_M-50_2016postVFP"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2016_Pt30M56to116"
        # eosDirSS        = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2016_Pt30M56to116_SS"
    elif args.beginRun > 297020 and args.beginRun < 306828:    # 2017
        byLsCSV         = resPath+"/FillByLs_2017.csv"
        mcFile          = resPath+"/mcCorrections/corrections_nPU_2017.p"
        sigTemplates    = resPath+"/sigTemplates/DYJetsToLL_M-50_2017"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2017_Pt30M56to116"
        # eosDirSS        = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2017_Pt30M56to116_SS"
    elif args.beginRun >= 306926 and args.beginRun < 307083: # 2017 H
        byLsCSV = resPath+"/FillByLs_2017_lowPU.csv"
        mcFile          = None
        # sigTemplates    = resPath+"/sigTemplates/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIISummer20UL17MiniAOD-FlatPU0to75_106X_mc2017_realistic_v6-v2"
        sigTemplates    = resPath+"/sigTemplates/DYJetsToLL_M-50_2017H"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2017H_Pt30"
        # eosDirSS        = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2017H_Pt30M56to116_SS"
    elif args.beginRun >= 315252:                           # 2018
        byLsCSV         = resPath+"/FillByLs_2018.csv"
        mcFile          = resPath+"/mcCorrections/corrections_nPU_2018.p"
        sigTemplates    = resPath+"/sigTemplates/DYJetsToLL_M-50_2018"
        eosDir          = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2018_Pt30M56to116"
        # eosDirSS        = "/nfs/dust/cms/user/dwalter/data/Lumi/V11/Histograms_UL2018_Pt30M56to116_SS"
    else:
        byLsCSV         = None
        mcFile          = None
        sigTemplates    = None
        eosDir          = None
        # eosDirSS        = None

    eosDir          = eosDir        if eosDir       != None and eosDir       != "None" else args.dirDQM
    # eosDirSS        = eosDirSS      if eosDirSS     != None and eosDirSS     != "None" else args.bkgSubtract
    byLsCSV         = byLsCSV       if byLsCSV      != None and byLsCSV      != "None" else args.byLsCSV
    sigTemplates    = sigTemplates  if sigTemplates != None and sigTemplates != "None" else args.sigTemplates
    mcFile          = mcFile        if mcFile       != None and mcFile       != "None" else args.mcCorrections

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
            'sigmod_hlt':[2,2,2,2],
            'sigtempl_hlt':sigTemplates+"/template_HLT.root",
            'sigmod_sel':[2,2,2,2],
            'sigtempl_sel':sigTemplates+"/template_Sel.root",
            'sigmod_trk':[2,2,2,2],
            'sigtempl_trk':sigTemplates+"/template_Trk.root",
            'sigmod_sta':[2,2,2,2],
            'sigtempl_sta':sigTemplates+"/template_Sta.root",
            }
        fit_options.update(_optns)

    if args.bkgTemplates and args.bkgTemplates != "None":

        _optns = {
            "bkgmod_hlt": [7,7,7,7],
            "bkgshape_hlt": args.bkgTemplates,
            "bkgmod_sel": [7,7,7,7],
            "bkgshape_sel": args.bkgTemplates,
            "bkgmod_trk": [7,7,7,7],
            "bkgshape_trk": args.bkgTemplates,
            "bkgmod_sta": [7,7,7,7],
            "bkgshape_sta": args.bkgTemplates,
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

    if mcFile:
        print("Use MC corrections from: ")
        print("file: "+mcFile)
        corr = getMCCorrection(mcFile)
    else:
        corr = None

    ########### Constant settings ##########
    secPerLS = float(23.3)
    currentYear = 2017

    maximumLS = 5000
    LSperMeasurement = 100  # required number of lumi sections per measurement
    LumiPerMeasurement = 20  # minimum recorded lumi for one measurement in pb-1

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
    MassBin_ = int(MassMax_ - MassMin_)

    ROOT.set_massRange(MassMin_, MassMax_)

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
        # for subtraction
        hYieldBB_neg = ROOT.TH1D("h_mass_yieldBB_Negative","",MassBin_, MassMin_, MassMax_)
        hYieldBE_neg = ROOT.TH1D("h_mass_yieldBE_Negative","",MassBin_, MassMin_, MassMax_)
        hYieldEE_neg = ROOT.TH1D("h_mass_yieldEE_Negative","",MassBin_, MassMin_, MassMax_)
        hYieldBB_pos = ROOT.TH1D("h_mass_yieldBB_Positive","",MassBin_, MassMin_, MassMax_)
        hYieldBE_pos = ROOT.TH1D("h_mass_yieldBE_Positive","",MassBin_, MassMin_, MassMax_)
        hYieldEE_pos = ROOT.TH1D("h_mass_yieldEE_Positive","",MassBin_, MassMin_, MassMax_)

        # For fitting
        hYieldBB = ROOT.TH1D("h_mass_yieldBB_Z","",MassBin_, MassMin_, MassMax_)
        hYieldBE = ROOT.TH1D("h_mass_yieldBE_Z","",MassBin_, MassMin_, MassMax_)
        hYieldEE = ROOT.TH1D("h_mass_yieldEE_Z","",MassBin_, MassMin_, MassMax_)
        hHLTpassB = ROOT.TH1D("h_mass_HLT_pass_central","",MassBin_, MassMin_, MassMax_)
        hHLTfailB = ROOT.TH1D("h_mass_HLT_fail_central","",MassBin_, MassMin_, MassMax_)
        hHLTpassE = ROOT.TH1D("h_mass_HLT_pass_forward","",MassBin_, MassMin_, MassMax_)
        hHLTfailE = ROOT.TH1D("h_mass_HLT_fail_forward","",MassBin_, MassMin_, MassMax_)
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

            hYieldBB_neg.Add(load("h_mass_yieldBB_Negative"))
            hYieldBE_neg.Add(load("h_mass_yieldBE_Negative"))
            hYieldEE_neg.Add(load("h_mass_yieldEE_Negative"))
            hYieldBB_pos.Add(load("h_mass_yieldBB_Positive"))
            hYieldBE_pos.Add(load("h_mass_yieldBE_Positive"))
            hYieldEE_pos.Add(load("h_mass_yieldEE_Positive"))

            hYieldBB.Add(load("h_mass_yieldBB_Z"))
            hYieldBE.Add(load("h_mass_yieldBE_Z"))
            hYieldEE.Add(load("h_mass_yieldEE_Z"))
            hHLTpassB.Add(load("h_mass_HLT_pass_central"))
            hHLTfailB.Add(load("h_mass_HLT_fail_central"))
            hHLTpassE.Add(load("h_mass_HLT_pass_forward"))
            hHLTfailE.Add(load("h_mass_HLT_fail_forward"))
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

        # if args.bkgSubtract:
        #     result = subtraction(
        #         hYieldBB, hYieldBE, hYieldEE,
        #         hHLTpassB,hHLTfailB,hHLTpassE,hHLTfailE,
        #         hSITpassB,hSITfailB,hSITpassE,hSITfailE,
        #         hTrkpassB,hTrkfailB,hTrkpassE,hTrkfailE,
        #         hStapassB,hStafailB,hStapassE,hStafailE
        #         )
        # else:
        result = measurement(0,
            hYieldBB, hYieldBE, hYieldEE,
            hYieldBB_neg, hYieldBE_neg, hYieldEE_neg,
            hYieldBB_pos, hYieldBE_pos, hYieldEE_pos,
            hHLTpassB, hHLTfailB, hHLTpassE, hHLTfailE,
            hSITpassB, hSITfailB, hSITpassE, hSITfailE,
            hTrkpassB, hTrkfailB, hTrkpassE, hTrkfailE,
            hStapassB, hStafailB, hStapassE, hStafailE,
            **fit_options
            )

    log.info("===Looping over runs...")
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
                if not mergeMeasurements and (recLumi >= LumiPerMeasurement or len(goodLSlist) >= LSperMeasurement):
                    break

            ### load histograms
            log.debug("Openning DQMIO.root file: %s", eosFile)

            def load(name_):
                return load_histo(name_, eosFile, goodLSlist, prefix_dqm, run, "new")

            if not args.inclusive:

                hYieldBB_neg = load("h_mass_yieldBB_Negative")
                hYieldBE_neg = load("h_mass_yieldBE_Negative")
                hYieldEE_neg = load("h_mass_yieldEE_Negative")
                hYieldBB_pos = load("h_mass_yieldBB_Positive")
                hYieldBE_pos = load("h_mass_yieldBE_Positive")
                hYieldEE_pos = load("h_mass_yieldEE_Positive")

                hYieldBB = load("h_mass_yieldBB_Z")
                hYieldBE = load("h_mass_yieldBE_Z")
                hYieldEE = load("h_mass_yieldEE_Z")
                hHLTpassB = load("h_mass_HLT_pass_central")
                hHLTfailB = load("h_mass_HLT_fail_central")
                hHLTpassE = load("h_mass_HLT_pass_forward")
                hHLTfailE = load("h_mass_HLT_fail_forward")
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

                # if args.bkgSubtract:
                #     # Do background subtraction from same sign events
                #     result = subtraction(
                #         hYieldBB, hYieldBE, hYieldEE,
                #         hHLTpassB,hHLTfailB,hHLTpassE,hHLTfailE,
                #         hSITpassB,hSITfailB,hSITpassE,hSITfailE,
                #         hTrkpassB,hTrkfailB,hTrkpassE,hTrkfailE,
                #         hStapassB,hStafailB,hStapassE,hStafailE
                #         )
                # else:
                # Do fits
                result = measurement(len(results),
                    hYieldBB, hYieldBE, hYieldEE,
                    hYieldBB_neg, hYieldBE_neg, hYieldEE_neg,
                    hYieldBB_pos, hYieldBE_pos, hYieldEE_pos,
                    hHLTpassB, hHLTfailB, hHLTpassE, hHLTfailE,
                    hSITpassB, hSITfailB, hSITpassE, hSITfailE,
                    hTrkpassB, hTrkfailB, hTrkpassE, hTrkfailE,
                    hStapassB, hStafailB, hStapassE, hStafailE,
                    **fit_options
                    )

            if result is None:
                continue

            df = data_run.loc[data_run['ls'].isin(goodLSlist)]

            dqmfile = ROOT.TFile(eosFile)

            hist_ZyieldBB = dqmfile.Get(prefix_dqm+"h_mass_yieldBB_Z")
            hist_ZyieldBE = dqmfile.Get(prefix_dqm+"h_mass_yieldBE_Z")
            hist_ZyieldEE = dqmfile.Get(prefix_dqm+"h_mass_yieldEE_Z")

            result.update(ls_corrections(df, result, len(results), hist_ZyieldBB, hist_ZyieldBE, hist_ZyieldEE))

            results.append(result)

        ## Write per measurement csv file - one per run
        log.info("===Writing per Run CSV file")
        results = pd.concat([pd.DataFrame([result]) for result in results])

        with open(outCSVDir + '/csvfile{0}.csv'.format(run), 'w') as file:
            results.to_csv(file, index=False)


    if args.writeSummaryCSV:
        writeSummaryCSV(outCSVDir)

    log.info("===Done")
