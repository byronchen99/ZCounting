import ROOT
import pandas as pd
import glob
import os
import numpy as np
import json
import pdb
import uncertainties as unc
import gc

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime, getMCCorrection, unorm, pquad

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None

# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

### Helper functions ###
# ---------------------->>>

# <<<------------------------
def get_correlationfactors(
    data,
    res="res/correlationfactors.json",
    byLs = False
):

    with open(res) as file:
        res = json.load(file)

    retvals = []
    for region in ("BB","BE","EE"):
        params = res[era.replace("2017H", "2017")][region]["nominal"]

        factor = pquad(data['avgpu'], *params)

        if not byLs:
            # sum factors weighted with the recorded luminosity in this lumisection
            factor = sum(factor * data['recorded(/pb)'].values) / sum(data['recorded(/pb)'].values)

        retvals.append(factor)

    return retvals


# <<<------------------------
def do_fits(m,
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

    ROOT.calculateGloEfficiency(h1TrkBPass, h1TrkBFail, h1StaBFail,
        m, 0, sigmod_trk[0], bkgmod_trk[0], sigmod_trk[1], bkgmod_trk[1],
        0, sigtempl_trk, bkgshape_trk)
    ROOT.calculateGloEfficiency(h1TrkEPass, h1TrkEFail, h1StaEFail,
        m, 1, sigmod_trk[2], bkgmod_trk[2], sigmod_trk[3], bkgmod_trk[3],
        0, sigtempl_trk, bkgshape_trk)

    ROOT.calculateDataEfficiency(h1SITBPass, h1SITBFail,
        m, "Sel", 0, sigmod_sel[0], bkgmod_sel[0], sigmod_sel[1], bkgmod_sel[1],
        0, sigtempl_sel, bkgshape_sel)
    ROOT.calculateDataEfficiency(h1SITEPass, h1SITEFail,
        m, "Sel", 1, sigmod_sel[2], bkgmod_sel[2], sigmod_sel[3], bkgmod_sel[3],
        0, sigtempl_sel, bkgshape_sel)

    # Fit BB, BE and EE separately but fit 1HLT and 2HLT each time together
    ROOT.calculateHLTEfficiencyAndYield(h2HLTBB_, h1HLTBB_,
        m, "BB", sigmod_yield, bkgmod_yield, sigmod_yield, bkgmod_yield,
        cBB_, 0, sigtempl_yield)
    ROOT.calculateHLTEfficiencyAndYield(h2HLTBE_, h1HLTBE_,
        m, "BE", sigmod_yield, bkgmod_yield, sigmod_yield, bkgmod_yield,
        cBE_, 0, sigtempl_yield)
    ROOT.calculateHLTEfficiencyAndYield(h2HLTEE_, h1HLTEE_,
        m, "EE", sigmod_yield, bkgmod_yield, sigmod_yield, bkgmod_yield,
        cEE_, 0, sigtempl_yield)

def extract_results(directory, m, cBB=1, cBE=1, cEE=1):
    print("INFO: ===Extracting fit results")

    fGloB = ROOT.TFile(directory+"/workspace_Glo_B_{0}.root".format(m),"READ")
    rGloB = fGloB.Get("fitResult")
    wGloB = fGloB.Get("workspace")

    effTrkB = wGloB.var("effTrk")
    effStaB = wGloB.var("effSta")
    (TrkeffB_m, StaeffB_m) = unc.correlated_values_norm(
        [(effTrkB.getVal(), effTrkB.getError()),(effStaB.getVal(), effStaB.getError())],
        [[1.,rGloB.correlation("effTrk","effSta")],[rGloB.correlation("effTrk","effSta"),1.]])
    GloeffB_m = TrkeffB_m*StaeffB_m

    fGloE = ROOT.TFile(directory+"/workspace_Glo_E_{0}.root".format(m),"READ")
    rGloE = fGloE.Get("fitResult")
    wGloE = fGloE.Get("workspace")

    effTrkE = wGloE.var("effTrk")
    effStaE = wGloE.var("effSta")
    (TrkeffE_m, StaeffE_m) = unc.correlated_values_norm(
        [(effTrkE.getVal(), effTrkE.getError()),(effStaE.getVal(), effStaE.getError())],
        [[1.,rGloE.correlation("effTrk","effSta")],[rGloE.correlation("effTrk","effSta"),1.]])
    GloeffE_m = TrkeffE_m*StaeffE_m

    fSelB = ROOT.TFile(directory+"/workspace_Sel_B_{0}.root".format(m),"READ")
    wSelB = fSelB.Get("workspace")
    SeleffB_m = unc.ufloat(wSelB.var("eff").getVal(), wSelB.var("eff").getError())

    fSelE = ROOT.TFile(directory+"/workspace_Sel_E_{0}.root".format(m),"READ")
    wSelE = fSelE.Get("workspace")
    SeleffE_m = unc.ufloat(wSelE.var("eff").getVal(), wSelE.var("eff").getError())

    fBB = ROOT.TFile(directory+"/workspace_BB_{0}.root".format(m),"READ")
    wBB = fBB.Get("workspace")
    HLTeffBB_m  = unc.ufloat(wBB.var("eff").getVal(), wBB.var("eff").getError())
    zRecBB_m = unc.ufloat(wBB.var("Nsig").getVal(), wBB.var("Nsig").getError())
    NbkgHLT1BB_m = unc.ufloat(wBB.var("NbkgFail").getVal(), wBB.var("NbkgFail").getError())
    NbkgHLT2BB_m = unc.ufloat(wBB.var("NbkgPass").getVal(), wBB.var("NbkgPass").getError())

    fBE = ROOT.TFile(directory+"/workspace_BE_{0}.root".format(m),"READ")
    wBE = fBE.Get("workspace")
    HLTeffBE_m = unc.ufloat(wBE.var("eff").getVal(), wBE.var("eff").getError())
    zRecBE_m = unc.ufloat(wBE.var("Nsig").getVal(), wBE.var("Nsig").getError())
    NbkgHLT1BE_m = unc.ufloat(wBE.var("NbkgFail").getVal(), wBE.var("NbkgFail").getError())
    NbkgHLT2BE_m = unc.ufloat(wBE.var("NbkgPass").getVal(), wBE.var("NbkgPass").getError())

    fEE = ROOT.TFile(directory+"/workspace_EE_{0}.root".format(m),"READ")
    wEE = fEE.Get("workspace")
    HLTeffEE_m  = unc.ufloat(wEE.var("eff").getVal(), wEE.var("eff").getError())
    zRecEE_m = unc.ufloat(wEE.var("Nsig").getVal(), wEE.var("Nsig").getError())
    NbkgHLT1EE_m = unc.ufloat(wEE.var("NbkgFail").getVal(), wEE.var("NbkgFail").getError())
    NbkgHLT2EE_m = unc.ufloat(wEE.var("NbkgPass").getVal(), wEE.var("NbkgPass").getError())

    # ZtoMuMu efficiency purely from data
    ZBBeff_m = StaeffB_m * StaeffB_m * GloeffB_m * GloeffB_m
    ZBEeff_m = StaeffB_m * StaeffE_m * GloeffB_m * GloeffE_m
    ZEEeff_m = StaeffE_m * StaeffE_m * GloeffE_m * GloeffE_m

    res = {
        "SeleffB": SeleffB_m,
        "SeleffE": SeleffE_m,
        "GloeffB": GloeffB_m,
        "GloeffE": GloeffE_m,
        "TrkeffB": TrkeffB_m,
        "TrkeffE": TrkeffE_m,
        "StaeffB": StaeffB_m,
        "StaeffE": StaeffE_m,
        "HLTeffBB": HLTeffBB_m,
        "HLTeffBE": HLTeffBE_m,
        "HLTeffEE": HLTeffEE_m,
        "zRecBB": zRecBB_m,
        "zRecBE": zRecBE_m,
        "zRecEE": zRecEE_m,
        "ZBBeff": ZBBeff_m,
        "ZBEeff": ZBEeff_m,
        "ZEEeff": ZEEeff_m,
        "zDelBB": zRecBB_m / ZBBeff_m,
        "zDelBE": zRecBE_m / ZBEeff_m,
        "zDelEE": zRecEE_m / ZEEeff_m,
        "NbkgHLT1BB": NbkgHLT1BB_m,
        "NbkgHLT1BE": NbkgHLT1BE_m,
        "NbkgHLT1EE": NbkgHLT1EE_m,
        "NbkgHLT2BB": NbkgHLT2BB_m,
        "NbkgHLT2BE": NbkgHLT2BE_m,
        "NbkgHLT2EE": NbkgHLT2EE_m,
        "cBB": cBB,
        "cBE": cBE,
        "cEE": cEE,
        "SeleffB_chi2pass": wSelB.arg("chi2pass").getVal(),
        "SeleffB_chi2fail": wSelB.arg("chi2fail").getVal(),
        "SeleffE_chi2pass": wSelE.arg("chi2pass").getVal(),
        "SeleffE_chi2fail": wSelE.arg("chi2fail").getVal(),
        "GloeffB_chi2pass": wGloB.arg("chi2pass").getVal(),
        "GloeffE_chi2pass": wGloE.arg("chi2pass").getVal(),
        "TrkeffB_chi2fail": wGloB.arg("chi2failTrk").getVal(),
        "TrkeffE_chi2fail": wGloE.arg("chi2failTrk").getVal(),
        "StaeffB_chi2fail": wGloB.arg("chi2failSta").getVal(),
        "StaeffE_chi2fail": wGloE.arg("chi2failSta").getVal(),
        "yieldBB_chi2pass": wBB.arg("chi2pass").getVal(),
        "yieldBB_chi2fail": wBB.arg("chi2fail").getVal(),
        "yieldBE_chi2pass": wBE.arg("chi2pass").getVal(),
        "yieldBE_chi2fail": wBE.arg("chi2fail").getVal(),
        "yieldEE_chi2pass": wEE.arg("chi2pass").getVal(),
        "yieldEE_chi2fail": wEE.arg("chi2fail").getVal(),
        "SeleffB_chi2": wSelB.arg("chi2").getVal(),
        "SeleffE_chi2": wSelE.arg("chi2").getVal(),
        "GloeffB_chi2": wGloB.arg("chi2").getVal(),
        "GloeffE_chi2": wGloE.arg("chi2").getVal(),
        "yieldBB_chi2": wBB.arg("chi2").getVal(),
        "yieldBE_chi2": wBE.arg("chi2").getVal(),
        "yieldEE_chi2": wEE.arg("chi2").getVal()
    }

    fGloB.Close()
    fGloB.Delete()
    wGloB.Delete()
    rGloB.Delete()

    fGloE.Close()
    fGloE.Delete()
    wGloE.Delete()
    rGloE.Delete()

    fSelB.Close()
    fSelB.Delete()
    wSelB.Delete()

    fSelE.Close()
    fSelE.Delete()
    wSelE.Delete()

    fBB.Close()
    fBB.Delete()
    wBB.Delete()

    fBE.Close()
    fBE.Delete()
    wBE.Delete()

    fEE.Close()
    fEE.Delete()
    wEE.Delete()

    return res

# <<<------------------------
def ls_corrections(data_m, result_m, m,
    h2HLTBB, h2HLTBE, h2HLTEE,
    h1HLTBB, h1HLTBE, h1HLTEE
):
    data_m["cBB"],data_m["cBE"],data_m["cEE"] = get_correlationfactors(data_m, correlationsFile, byLs=True)

    # --- per ls dataframe, one per measurement
    for eta, h1, h2 in (("BB",h1HLTBB, h2HLTBB), ("BE",h1HLTBE, h2HLTBE), ("EE",h1HLTEE, h2HLTEE)):

        zRec_m = result_m['zRec{}'.format(eta)]
        corr_m = result_m['c{}'.format(eta)]
        HLTeff_m = result_m['HLTeff{}'.format(eta)]

        NsigHLT1_m = 2 * zRec_m * HLTeff_m * (1 - HLTeff_m * corr_m)
        NsigHLT2_m = zRec_m * HLTeff_m**2 * corr_m

        NbkgHLT1_m = result_m['NbkgHLT1{}'.format(eta)]
        NbkgHLT2_m = result_m['NbkgHLT2{}'.format(eta)]

        purity1_m = NsigHLT1_m / (NsigHLT1_m+NbkgHLT1_m)
        purity2_m = NsigHLT2_m / (NsigHLT2_m+NbkgHLT2_m)

        zEff_m = result_m["Z{}eff".format(eta)]

        data_m['N1HLT{}'.format(eta)] = data_m['ls'].apply(lambda ls, h=h1: h.ProjectionY("", ls, ls, "e").Integral(MassMin_-30, MassMax_-30))
        data_m['N2HLT{}'.format(eta)] = data_m['ls'].apply(lambda ls, h=h2: h.ProjectionY("", ls, ls, "e").Integral(MassMin_-30, MassMax_-30))

        data_m['Nsig1HLT{}'.format(eta)] = purity1_m.n * data_m['N1HLT{}'.format(eta)]
        data_m['Nsig2HLT{}'.format(eta)] = purity2_m.n * data_m['N2HLT{}'.format(eta)]

        # eff_HLT = 2*n2 /(c*(n1 + 2*n2))
        data_m['HLTeff{}_mc'.format(eta)] = 2*data_m['Nsig2HLT{}'.format(eta)] / (data_m['c{}'.format(eta)] * (data_m['Nsig1HLT{}'.format(eta)] + 2*data_m['Nsig2HLT{}'.format(eta)]))

        # nZ = (n2 + n1/2.)**2 * c / n2
        data_m['zRec{}_mc'.format(eta)] = (data_m['Nsig2HLT{}'.format(eta)] + data_m['Nsig1HLT{}'.format(eta)]/2.)**2 / data_m['Nsig2HLT{}'.format(eta)]

        df.loc[~np.isfinite(df['zRec{}_mc'.format(eta)]), 'zRec{}_mc'.format(eta)] = 0

        data_m['zRec{}_mcUp'.format(eta)]   = data_m['zRec{}_mc'.format(eta)] * (2 * (data_m['c{}'.format(eta)]-1) + 1)
        data_m['zRec{}_mcDown'.format(eta)] = data_m['zRec{}_mc'.format(eta)] * (0.5 * (data_m['c{}'.format(eta)]-1) + 1)
        data_m['zRec{}_mc'.format(eta)]     = data_m['zRec{}_mc'.format(eta)] * data_m['c{}'.format(eta)]

        data_m['zDel{}_mc'.format(eta)] = data_m['zRec{}_mc'.format(eta)] / zEff_m.n
        data_m['zDel{}_mcUp'.format(eta)] = data_m['zRec{}_mcUp'.format(eta)] / zEff_m.n
        data_m['zDel{}_mcDown'.format(eta)] = data_m['zRec{}_mcDown'.format(eta)] / zEff_m.n

    # sort out lumisections without any Z candidate (maybe trigger was off)
    # consider lumisections where we would expect to have at least any count (> 0.01 /pb we expect 0.01*500 = 5 Z bosons)
    zeros = (data_m['N1HLTBB']+data_m['N2HLTBB']+data_m['N1HLTBE']+data_m['N2HLTBE']+data_m['N1HLTEE']+data_m['N2HLTEE'] == 0) & (data_m['recorded(/pb)'] > 0.01)
    if sum(zeros) > 0:
        print("sort out {0} lumisections in run {1}: {2}".format(sum(zeros), run, data_m[zeros]['ls'].values))
        data_m = data_m[~zeros]

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
        "zRecBB_mc": data_m['zRecBB_mc'].sum(),
        "zRecBE_mc": data_m['zRecBE_mc'].sum(),
        "zRecEE_mc": data_m['zRecEE_mc'].sum(),
        "zDelBB_mc": data_m['zDelBB_mc'].sum(),
        "zDelBE_mc": data_m['zDelBE_mc'].sum(),
        "zDelEE_mc": data_m['zDelEE_mc'].sum(),
        "zDelBB_mcUp": data_m['zDelBB_mcUp'].sum(),
        "zDelBE_mcUp": data_m['zDelBE_mcUp'].sum(),
        "zDelEE_mcUp": data_m['zDelEE_mcUp'].sum(),
        "zDelBB_mcDown": data_m['zDelBB_mcDown'].sum(),
        "zDelBE_mcDown": data_m['zDelBE_mcDown'].sum(),
        "zDelEE_mcDown": data_m['zDelEE_mcDown'].sum(),
        "lumiDel": delLumi_m,
        "lumiRec": recLumi_m,
        "timewindow": timeWindow_m,
        "deadtime": deadtime_m,
        "pileUp": data_m['avgpu'].mean(),
        # "ZBBeff": data_m['effBB'].values.mean(),
        # "ZBEeff": data_m['effBE'].values.mean(),
        # "ZEEeff": data_m['effEE'].values.mean(),
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

    print("INFO: ===Writing overall CSV file")
    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile??????.csv'))
    df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)

    with open(outCSVDir + '/Mergedcsvfile_perMeasurement.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

    print("INFO: ===Writing overall CSV file per LS")
    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile*_*.csv'))
    csvList = []
    # add measurement label to the csv list
    for m in rateFileList:
        csv = pd.read_csv(m)
        measurement = int(m.split("_")[-1][:-4])
        fill = int(m.split("_")[-2].split("csvfile")[-1])
        csv["measurement"] = measurement
        csvList.append(csv)
    df_merged = pd.concat(csvList, ignore_index=True)

    # df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)

    with open(outCSVDir + '/Mergedcsvfile_perLS.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

## Get root file with histograms
def getFileName(directory, run):
    # check if run was processed already
    print("INFO: ===Checking input DQMIO.root file...")
    eosFileList = glob.glob(directory + '/*/*' + str(run) + '*root')
    if not len(eosFileList) > 0:
        print("INFO: The file does not yet exist for run: " + str(run))
        return None
    elif len(eosFileList) > 1:
        print("INFO: Multiple files found for run: " + str(run))
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
    parser.add_argument("-m", "--measurement", help="Only fit a specific measurement of the run", type=int, default=None)
    parser.add_argument('--mcCorrections', default="default", type=str,
                        help='specify .json file with MC corrections for muon correlations')
    parser.add_argument("-c", "--writeSummaryCSV", default=False, action="store_true",
                        help="produce merged CSV with all runs")
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
    parser.add_argument('--collect', default=False, action="store_true",
                        help='specify whether or not to run the fits or just collect the results')
    parser.add_argument("-o", "--dirOut", help="where to store the output files", default="./")

    args = parser.parse_args()


    ########################################
    # link to resouces
    prefix_dqm="" #"DQMData/Run {0}/ZCounting/Run summary/Histograms/".format(run)
    resPath = cmsswbase + "/src/ZCounting/ZHarvester/res/"
    if( args.beginRun >= 272007 and args.beginRun < 278808
        # there is an overlap for 2016 F in runs with pre and post VFP settings
        and args.beginRun not in [278769, 278801, 278802, 278803, 278804, 278805, 278808]
    ): # 2016 pre VFP
        byLsCSV         = resPath+"/FillByLs_RunII.csv"
        correlationsFile          = resPath+"/correlations/CorrFactor-V14_04_02_quad-d20210909-t104351/info.json"
        sigTemplates    = resPath+"/sigTemplates/V13_02/ZSignalTemplate-V13_02-Summer16preVFP-DYJetsToLL_M_50_NLO.root"
        bkgTemplates    = resPath+"/qcdTemplates/Run2016preVFP.root"
        era = "2016preVFP"
    elif args.beginRun < 294645:    # 2016 post VFP
        byLsCSV         = resPath+"/FillByLs_RunII.csv"
        correlationsFile          = resPath+"/correlations/CorrFactor-V14_04_02_quad-d20210909-t104351/info.json"
        sigTemplates    = resPath+"/sigTemplates/V13_02/ZSignalTemplate-V13_02-Summer16postVFP-DYJetsToLL_M_50_NLO.root"
        bkgTemplates    = resPath+"/qcdTemplates/Run2016postVFP.root"
        era = "2016postVFP"
    elif args.beginRun > 297020 and args.beginRun < 306828:     # 2017
        byLsCSV         = resPath+"/FillByLs_RunII.csv"
        correlationsFile          = resPath+"/correlations/CorrFactor-V14_04_02_quad-d20210909-t104351/info.json"
        sigTemplates    = resPath+"/sigTemplates/V13_02/ZSignalTemplate-V13_02-Fall17-DYJetsToLL_M_50_NLO.root"
        bkgTemplates    = resPath+"/qcdTemplates/Run2017.root"
        era = "2017"
    elif args.beginRun >= 306926 and args.beginRun < 307083:    # 2017 H
        byLsCSV         = resPath+"/FillByLs_2017_lowPU.csv"
        correlationsFile          = resPath+"/correlations/CorrFactor-V14_04_02_SingleMu25_quad-d20210909-t104733/info.json"
        sigTemplates    = resPath+"/sigTemplates/V13_02/ZSignalTemplate-V13_02_SingleMu25-Fall17-DYJetsToLL_M_50_NLO.root"
        era = "2017H"
    elif args.beginRun >= 315252:                               # 2018
        byLsCSV         = resPath+"/FillByLs_RunII.csv"
        correlationsFile          = resPath+"/correlations/CorrFactor-V14_04_02_quad-d20210909-t104351/info.json"
        sigTemplates    = resPath+"/sigTemplates/V13_02/ZSignalTemplate-V13_02-Autumn18-DYJetsToLL_M_50_NLO.root"
        bkgTemplates    = resPath+"/qcdTemplates/Run2018.root"
        era = "2018"
    else:
        byLsCSV             = None
        correlationsFile    = None
        sigTemplates        = None

    eosDir           = args.dirDQM
    byLsCSV          = byLsCSV          if args.byLsCSV       == "default"   else args.byLsCSV
    sigTemplates     = sigTemplates     if args.sigTemplates  == "default"   else args.sigTemplates
    bkgTemplates     = bkgTemplates     if args.bkgTemplates  == "default"   else args.bkgTemplates
    correlationsFile = correlationsFile if args.mcCorrections == "default"   else args.mcCorrections
    measurement       = args.measurement

    print("----------------------------------")
    print("Use eosDir:              {0}".format(eosDir))
    print("Use byLsCSV:             {0}".format(byLsCSV))
    print("Use sigTemplates:        {0}".format(sigTemplates))
    print("Use bkgTemplates:        {0}".format(bkgTemplates))
    print("Use correlationsFile:    {0}".format(correlationsFile))
    print("Mass range from:         {0} to {1}".format(*args.mass))
    print("Lumi per Measurement:    {0}".format(args.LumiPerMeasurement))
    print("----------------------------------")

    ########## Input configuration ##########
    # ByLS csv inputs generated by testBRIL.sh
    byLS_filelist = glob.glob(byLsCSV)
    byLS_filelist.sort(key=os.path.getmtime)
    byLS_filename = byLS_filelist[-1]
    print("INFO: The brilcalc csv file: " + str(byLS_filename))

    # inputs: to build MC*Gaussian template for efficiency fitting
    fit_options={}
    if sigTemplates and sigTemplates != "None":
        _optns = {
            'sigmod_yield':2,
            'sigtempl_yield':sigTemplates,
            # 'sigmod_hlt':[2,2,2,2],
            # 'sigtempl_hlt':sigTemplates,
            'sigmod_sel':[2,2,2,2],
            'sigtempl_sel':sigTemplates,
            'sigmod_trk':[2,2,2,2],
            'sigtempl_trk':sigTemplates,
            'sigmod_sta':[2,2,2,2],
            'sigtempl_sta':sigTemplates,
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
            print("WARNING: directory already exists ...")

    # if correlationsFile:
    #     print("Use MC corrections from: ")
    #     print("file: "+correlationsFile)
    #     corr = getMCCorrection(correlationsFile)
    # else:
    #     corr = None

    ########### Constant settings ##########
    secPerLS = float(23.3)
    currentYear = 2017

    maximumLS = 5000
    LumiPerMeasurement = args.LumiPerMeasurement  # minimum recorded lumi for one measurement in pb-1

    #configuration of fit models
    MassMin_ = int(args.mass[0])
    MassMax_ = int(args.mass[1])
    MassBin_ = int(args.mass[2])

    if not args.collect:
        print("INFO: Loading C marco...")
        # load functions for fitting
        ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(
            __file__)) + "/calculateDataEfficiency.C")

        ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

        ROOT.set_ptCut(args.ptCut)

    print("INFO: Loading input byls csv...")
    byLS_file = open(str(byLS_filename))
    byLS_lines = byLS_file.readlines()
    byLS_data = pd.read_csv(byLS_filename, sep=',', low_memory=False,
        skiprows=lambda x: byLS_lines[x].startswith('#') and not byLS_lines[x].startswith('#run'))

    print("INFO: formatting csv file...")    # formatting the csv
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

    if not args.collect and args.inclusive:
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
                print("INFO: Continue")
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

        # get correlation factor for second muon
        cBB, cBE, cEE = get_correlationfactors(byLS_data, correlationsFile)

        do_fits(0,
            h2HLTBB, h2HLTBE, h2HLTEE,
            h1HLTBB, h1HLTBE, h1HLTEE,
            cBB, cBE, cEE,
            hSITpassB, hSITfailB, hSITpassE, hSITfailE,
            hTrkpassB, hTrkfailB, hTrkpassE, hTrkfailE,
            hStapassB, hStafailB, hStapassE, hStafailE,
            **fit_options
            )

    print("INFO: ===Looping over runs... {0} to {1}".format(int(args.beginRun), int(args.endRun)))
    for run, data_run in byLS_data.groupby('run'):
        if run < int(args.beginRun) or run >= int(args.endRun):
            continue

        # clean garbage collector (objects that were deleted with 'del')
        gc.collect()

        fill = data_run.drop_duplicates('fill')['fill'].values[0]
        LSlist = data_run.query('ls <= {0}'.format(maximumLS))['ls'].values.tolist()

        # Consider only runs with at least 5 LS
        if len(LSlist) <=5:
            continue

        outSubDir = outDir + "Run{0}/".format(run)
        if not os.path.isdir(outSubDir):
            os.mkdir(outSubDir)

        print("INFO: ===Running Fill {0}".format(fill))
        print("INFO: ===Running Run {0}".format(run))

        eosFile = getFileName(eosDir, run)      # >>> histograms binned in mass
        if eosFile is None:
            print("INFO: Continue")
            continue

        print("INFO: ===Looping over measurements...")
        results = []
        m = 0
        while len(LSlist) > 0:  # begin next measurement "m"

            # merge data to one measuement if remaining luminosity is too less for two measuements
            mergeMeasurements = sum(
                data_run.loc[data_run['ls'].isin(LSlist)]['recorded(/pb)'].values) < 1.5 * LumiPerMeasurement

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

            if measurement is not None and measurement != m:
                m += 1
                continue

            # get correlation factor for second muon
            cBB, cBE, cEE = get_correlationfactors(df, correlationsFile)

            if not args.collect and not args.inclusive:

                ### load histograms
                def load(name_):
                    return load_histo(name_, eosFile, goodLSlist, prefix_dqm, run, "new")

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

                do_fits(m,
                    h2HLTBB, h2HLTBE, h2HLTEE,
                    h1HLTBB, h1HLTBE, h1HLTEE,
                    cBB, cBE, cEE,
                    hSITpassB, hSITfailB, hSITpassE, hSITfailE,
                    hTrkpassB, hTrkfailB, hTrkpassE, hTrkfailE,
                    hStapassB, hStafailB, hStapassE, hStafailE,
                    **fit_options
                    )
            if args.inclusive:
                result = extract_results(outDir+"Fits_Inclusive/", 0, cBB, cBE, cEE)
            else:
                result = extract_results(outSubDir, m, cBB, cBE, cEE)

            if result:
                file_ = ROOT.TFile(eosFile)
                h2HLTBB_ls = file_.Get("h_mass_2HLTBB_Z")
                h1HLTBB_ls = file_.Get("h_mass_1HLTBB_Z")
                h2HLTBE_ls = file_.Get("h_mass_2HLTBE_Z")
                h1HLTBE_ls = file_.Get("h_mass_1HLTBE_Z")
                h2HLTEE_ls = file_.Get("h_mass_2HLTEE_Z")
                h1HLTEE_ls = file_.Get("h_mass_1HLTEE_Z")

                result.update(ls_corrections(df, result, m,
                    h2HLTBB_ls, h2HLTBE_ls, h2HLTEE_ls, h1HLTBB_ls, h1HLTBE_ls, h1HLTEE_ls))
                results.append(result)
            else:
                print("INFO: === No result - continue")

            m += 1

        ## Write per measurement csv file - one per run
        print("INFO: ===Writing per Run CSV file")
        results = pd.concat([pd.DataFrame([result]) for result in results])

        with open(outCSVDir + '/csvfile{0}.csv'.format(run), 'w') as file:
            results.to_csv(file, index=False)


    if args.writeSummaryCSV:
        writeSummaryCSV(outCSVDir)

    print("INFO: ===Done")
