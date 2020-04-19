

def full_efficiency(effBB, effBE, effEE):
    return effBB * ZBBAcc + effBE * ZBEAcc + effEE * ZEEAcc

def full_uncertainty(effBB_stat, effBE_stat, effEE_stat):
    return [np.sqrt((ZBBAcc * effBB_stat[i]) ** 2 + (ZBEAcc * effBE_stat[i]) ** 2 + (ZEEAcc * effEE_stat[i]) ** 2) for i in range(0,1)]

def measurement(recLumi, m, outputDir, h1Reco,
    h1HLTBPass, h1HLTBFail, h1HLTEPass, h1HLTEFail,
    h1SITBPass, h1SITBFail, h1SITEPass, h1SITEFail,
    h1GloBPass, h1GloBFail, h1GloEPass, h1GloEFail,
    sigmod_yield=1, bkgmod_yield=5,
    sigtempl_yield="", bkghist_yield=0,
    sigmod_hlt=[1,1,1,1], bkgmod_hlt=[5,5,5,5],
    sigtempl_hlt="", bkgshape_hlt="",
    sigmod_sel=[1,1,1,1], bkgmod_sel=[5,5,5,5],
    sigtempl_sel="", bkgshape_sel="",
    sigmod_glo=[1,1,1,1], bkgmod_glo=[5,5,5,5],
    sigtempl_glo="", bkgshape_glo="",
    ):
    ## compute Z yield #TODO
    Zyieldres_m = ROOT.getZyield(h1Reco, outputDir, m, sigmod_yield, bkgmod_yield,
        ptCutTag, ptCutProbe, recLumi, sigtempl_yield, bkghist_yield)


    ### compute muon efficiencies
    HLTeffresB_m = ROOT.calculateDataEfficiency(h1HLTBPass, h1HLTBFail,
                                                outputDir, m, "HLT", 0, sigmod_hlt[0], bkgmod_hlt[0], sigmod_hlt[1], bkgmod_hlt[1],
                                                ptCutTag, ptCutProbe, 0, recLumi, sigtempl_hlt, bkgshape_hlt)
    HLTeffresE_m = ROOT.calculateDataEfficiency(h1HLTEPass, h1HLTEFail,
                                                outputDir, m, "HLT", 1, sigmod_hlt[2], bkgmod_hlt[2], sigmod_hlt[3], bkgmod_hlt[3],
                                                ptCutTag, ptCutProbe, 0, recLumi, sigtempl_hlt, bkgshape_hlt)

    SeleffresB_m = ROOT.calculateDataEfficiency(h1SITBPass, h1SITBFail,
                                                outputDir, m, "Sel", 0, sigmod_sel[0], bkgmod_sel[0], sigmod_sel[1], bkgmod_sel[1],
                                                ptCutTag, ptCutProbe, 0, recLumi, sigtempl_sel, bkgshape_sel)
    SeleffresE_m = ROOT.calculateDataEfficiency(h1SITEPass, h1SITEFail,
                                                outputDir, m, "Sel", 1, sigmod_sel[2], bkgmod_sel[2], sigmod_sel[3], bkgmod_sel[3],
                                                ptCutTag, ptCutProbe, 0, recLumi, sigtempl_sel, bkgshape_sel)

    GloeffresB_m = ROOT.calculateDataEfficiency(h1GloBPass, h1GloBFail,
                                                outputDir, m, "Glo", 0, sigmod_glo[0], bkgmod_glo[0], sigmod_glo[1], bkgmod_glo[1],
                                                ptCutTag, ptCutProbe, 0, recLumi, sigtempl_glo, bkgshape_glo)
    GloeffresE_m = ROOT.calculateDataEfficiency(h1GloEPass, h1GloEFail,
                                                outputDir, m, "Glo", 1, sigmod_glo[2], bkgmod_glo[2], sigmod_glo[3], bkgmod_glo[3],
                                                ptCutTag, ptCutProbe, 0, recLumi, sigtempl_glo, bkgshape_glo)


    HLTeffB_m = HLTeffresB_m[0]
    HLTeffE_m = HLTeffresE_m[0]
    SeleffB_m = SeleffresB_m[0]
    SeleffE_m = SeleffresE_m[0]
    GloeffB_m = GloeffresB_m[0]
    GloeffE_m = GloeffresE_m[0]

    # ZtoMuMu efficiency purely from data
    ZBBEff = (GloeffB_m * GloeffB_m * SeleffB_m * SeleffB_m * (1 - (1 - HLTeffB_m) * (1 - HLTeffB_m)))
    ZBEEff = (GloeffB_m * GloeffE_m * SeleffB_m * SeleffE_m * (1 - (1 - HLTeffB_m) * (1 - HLTeffE_m)))
    ZEEEff = (GloeffE_m * GloeffE_m * SeleffE_m * SeleffE_m * (1 - (1 - HLTeffE_m) * (1 - HLTeffE_m)))

    # Statistic Uncertainties (low,high) error propagation
    ZBBEff_EStat = [0., 0.]
    ZBEEff_EStat = [0., 0.]
    ZEEEff_EStat = [0., 0.]
    for i in (1, 2):
        ZBBEff_EStat[i - 1] = 2 * ZBBEff * np.sqrt(
            (GloeffresB_m[i] / GloeffB_m) ** 2 +
            (SeleffresB_m[i] / SeleffB_m) ** 2 +
            ((1 - HLTeffB_m) / (1 - (1 - HLTeffB_m) ** 2) * HLTeffresB_m[i]) ** 2
        )
        ZEEEff_EStat[i - 1] = 2 * ZEEEff * np.sqrt(
            (GloeffresE_m[i] / GloeffE_m) ** 2 +
            (SeleffresE_m[i] / SeleffE_m) ** 2 +
            ((1 - HLTeffE_m) / (1 - (1 - HLTeffE_m) ** 2) * HLTeffresE_m[i]) ** 2
        )
        ZBEEff_EStat[i - 1] = ZBEEff * np.sqrt(
            (GloeffresB_m[i] / GloeffB_m) ** 2 +
            (GloeffresE_m[i] / GloeffE_m) ** 2 +
            (SeleffresB_m[i] / SeleffB_m) ** 2 +
            (SeleffresE_m[i] / SeleffE_m) ** 2 +
            ((1 - HLTeffE_m) / (1 - (1 - HLTeffB_m) * (1 - HLTeffE_m)) * HLTeffresB_m[i]) ** 2 +
            ((1 - HLTeffB_m) / (1 - (1 - HLTeffB_m) * (1 - HLTeffE_m)) * HLTeffresE_m[i]) ** 2
        )

    ZEff = full_efficiency(ZBBEff, ZBEEff, ZEEEff)
    ZEff_stat = full_uncertainty(ZBBEff_EStat, ZBEEff_EStat, ZEEEff_EStat)

    res = {
        "zYield": Zyieldres_m[0],
        "zYield_err": max(Zyieldres_m[1:3]),
        "zYield_chi2": Zyieldres_m[3],
        "zYield_purity": Zyieldres_m[4],
        "zYield_purity_err": max(Zyieldres_m[5:]),
        "HLTeffB": HLTeffB_m,
        "HLTeffE": HLTeffE_m,
        "SeleffB": SeleffB_m,
        "SeleffE": SeleffE_m,
        "GloeffB": GloeffB_m,
        "GloeffE": GloeffE_m,
        "HLTeffB_chi2pass": HLTeffresB_m[3],
        "HLTeffB_chi2fail": HLTeffresB_m[4],
        "HLTeffE_chi2pass": HLTeffresE_m[3],
        "HLTeffE_chi2fail": HLTeffresE_m[4],
        "SeleffB_chi2pass": SeleffresB_m[3],
        "SeleffB_chi2fail": SeleffresB_m[4],
        "SeleffE_chi2pass": SeleffresE_m[3],
        "SeleffE_chi2fail": SeleffresE_m[4],
        "GloeffB_chi2pass": GloeffresB_m[3],
        "GloeffB_chi2fail": GloeffresB_m[4],
        "GloeffE_chi2pass": GloeffresE_m[3],
        "GloeffE_chi2fail": GloeffresE_m[4],
        "Zeff": ZEff,
        "ZBBeff": ZBBEff,
        "ZBEeff": ZBEEff,
        "ZEEeff": ZEEEff,
        "Zeff_stat": max(ZEff_stat),
    }
    return res

def ls_corrections(data_m, result_m, m, h2Zyield, h1ZyieldBB, h1ZyieldEE):
    # --- per ls dataframe, one per measurement
    data_m['effBB_mc'] = result_m["ZBBeff"] - (data_m['avgpu'] * corr['BB_a'] + corr['BB_b'])
    data_m['effBE_mc'] = result_m["ZBEeff"] - (data_m['avgpu'] * corr['BE_a'] + corr['BE_b'])
    data_m['effEE_mc'] = result_m["ZEEeff"] - (data_m['avgpu'] * corr['EE_a'] + corr['EE_b'])
    data_m['eff_mc'] = full_efficiency(data_m['effBB_mc'], data_m['effBE_mc'], data_m['effEE_mc'])

    data_m['yield'] = data_m['ls'].apply(lambda ls: h2Zyield.ProjectionY("", ls, ls, "e").GetEntries())
    data_m['yieldBB'] = data_m['ls'].apply(lambda ls: h1ZyieldBB.GetBinContent(ls))
    data_m['yieldEE'] = data_m['ls'].apply(lambda ls: h1ZyieldEE.GetBinContent(ls))

    data_m['zYield'] = data_m['yield'] * result_m["zYield_purity"]
    data_m['zYieldBB'] = data_m['yieldBB'] * result_m["zYield_purity"]
    data_m['zYieldEE'] = data_m['yieldEE'] * result_m["zYield_purity"]

    data_m['zDel'] = data_m['zYield'] / result_m["Zeff"]
    data_m['zDelBB'] = data_m['zYieldBB'] / result_m["ZBBeff"]
    data_m['zDelEE'] = data_m['zYieldEE'] / result_m["ZEEeff"]

    data_m['zDel_mc'] = data_m['zYield'] / data_m['eff_mc']
    data_m['zDelBB_mc'] = data_m['zYieldBB'] / (data_m['effBB_mc'])
    data_m['zDelEE_mc'] = data_m['zYieldEE'] / (data_m['effEE_mc'])

    data_m['z_relstat'] = (1./np.sqrt(data_m['zYield']) \
        + result_m["Zeff_stat"]/result_m["Zeff"] \
        + result_m["zYield_purity_err"]/result_m["zYield_purity"]).replace(np.inf, 1.)

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

    zDel_m = data_m['zDel'].sum()
    zDel_mc_m = data_m['zDel_mc'].sum()
    z_relstat_m = result_m['zYield_err']/result_m['zYield'] + result_m["Zeff_stat"]/result_m["Zeff"]
    zRate_m = zDel_m / (timeWindow_m * deadtime_m)
    zRate_mc_m = zDel_mc_m / (timeWindow_m * deadtime_m)
    zLumi_m = zDel_m / sigma_fid
    zLumi_mc_m = zDel_mc_m / sigma_fid
    zXSec_m = zDel_m / recLumi_m
    zXSec_mc_m = zDel_mc_m / recLumi_m

    res = {
        "fill": fill,
        "run": run,
        "tdate_begin": min(data_m['time']),
        "tdate_end": max(data_m['time']),
        "yield": data_m['yield'],
        "zYield": data_m['zYield'],
        "zDel": zDel_m,
        "zDel_mc": zDel_mc_m,
        "z_relstat": z_relstat_m,
        "zRate": zRate_m,
        "zRate_mc": zRate_mc_m,
        "zLumi": zLumi_m,
        "zLumi_mc": zLumi_mc_m,
        "zXSec": zXSec_m,
        "zXSec_mc": zXSec_mc_m,
        "lumiDel": delLumi_m,
        "lumiRec": recLumi_m,
        "timewindow": timeWindow_m,
        "deadtime": deadtime_m,
        "pileUp": data_m['avgpu'].mean(),
        "ZMCeff": data_m['eff_mc'].mean(),
        "ZMCeffBB": data_m['effBB_mc'].mean(),
        "ZMCeffBE": data_m['effBE_mc'].mean(),
        "ZMCeffEE": data_m['effEE_mc'].mean(),
    }
    return res

def load_histo(name_, dqmfile_, lumisections_=[0,], prefix_="", run_=0):

    h_X_ls = dqmfile_.Get("{0}{1}".format(prefix_,name_))
    h_X = h_X_ls.ProjectionY("h_mass_{0}_{1}".format(name_, run_), lumisections_[0], lumisections_[0], "e")
    for ls in lumisections_[1:]:
        h_X.Add(h_X_ls.ProjectionY("h_mass_{0}_{1}_{2}".format(name_, run_, ls), ls, ls, "e"))
    return h_X


## Write big per measurement csv file
def writeSummaryCSV(outCSVDir):
    log.info("===Writing overall CSV file")
    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile??????.csv'))
    df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)
    log.info("total number of delivered Z is: {0} fb".format(sum(df_merged['zDel_mc'])))
    log.info("total fiducial cross section is: {0} fb".format(sum(df_merged['zDel_mc'])/sum(df_merged['lumiRec'])))

    with open(outCSVDir + 'Mergedcsvfile_perMeasurement.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile*_*.csv'))
    df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)

    with open(outCSVDir + 'Mergedcsvfile_perLS.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

################################################################################
if __name__ == '__main__':
    import argparse
    import logging as log
    import ROOT
    import pandas as pd
    import glob
    import os
    import numpy as np
    import json
    import pdb

    from Utils.Utils import to_RootTime

    # disable panda warnings when assigning a new column in the dataframe
    pd.options.mode.chained_assignment = None
    # turn off graphical output on screen
    ROOT.gROOT.SetBatch(True)
    cmsswbase = os.environ['CMSSW_BASE']

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beginRun", help="first run to analyze [%default]", default=272007)
    parser.add_argument("-e", "--endRun", help="analyze stops when comes to this run [%default]", default=1000000)
    parser.add_argument('--mcCorrections', default='./Resources/MCCorrections.json', type=str,
                        help='specify .json file with MC corrections for muon correlations')
    parser.add_argument("-v", "--verbose", help="increase logging level from INFO to DEBUG", default=False,
                        action="store_true")
    parser.add_argument("-c", "--writeSummaryCSV", help="produce merged CSV with all runs", default=True)
    parser.add_argument("-i", "--dirDQM", help="Directory to the input root files from the DQM Offline module",
                        default="/eos/home-d/dwalter/www/ZCounting/DQM-Offline-2018/")
    parser.add_argument("--byLsCSV", help="ByLs csv input generated by testBril.sh",
                        default="/nfs/dust/cms/user/dwalter/data/Lumi/For_2017/FillByLs_2017.csv")
    parser.add_argument("--sigTemplates", help="Directory to root file for Z mass template",
                        default=None, type=str)
    parser.add_argument("--bkgTemplates", help="Directory to root file for bkg template",
                        default=None, type=str)
    parser.add_argument('--ptCut', type=float, default=30.,
                        help='specify lower pt cut on tag and probe muons')
    parser.add_argument('--inclusive', default=False, action="store_true",
                        help='specify whether or not to do an inclusive fit of the specified runs')
    parser.add_argument("-o", "--dirOut", help="where to store the output files", default="./")

    args = parser.parse_args()
    if args.verbose:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
    else:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.INFO)

    ########## Input configuration ##########
    # ByLS csv inputs generated by testBRIL.sh
    byLS_filelist = glob.glob(args.byLsCSV)
    byLS_filelist.sort(key=os.path.getmtime)
    byLS_filename = byLS_filelist[-1]
    log.info("The brilcalc csv file: " + str(byLS_filename))

    eosDir = args.dirDQM
    prefix_dqm="" #"DQMData/Run {0}/ZCounting/Run summary/Histograms/".format(run)

    # inputs: to build MC*Gaussian template for efficiency fitting
    fit_options={}
    if args.sigTemplates and args.sigTemplates != "None":
        _optns = {
            'sigmod_yield':2,
            'sigtempl_yield':args.sigTemplates+"/template_ZYield.root",
            'sigmod_hlt':[2,2,2,2],
            'sigtempl_hlt':args.sigTemplates+"/template_HLT.root",
            'sigmod_sel':[2,2,2,2],
            'sigtempl_sel':args.sigTemplates+"/template_Sel.root",
            'sigmod_glo':[2,2,2,2],
            'sigtempl_glo':args.sigTemplates+"/template_Glo.root",
            }
        fit_options.update(_optns)

    if args.bkgTemplates and args.bkgTemplates != "None":
        tfile = ROOT.TFile.Open(args.bkgTemplates+"/Reco.root","READ")
        hBkg_yield = tfile.Get("KDE_cosine")
        hBkg_yield.SetDirectory(0)
        tfile.Close()

        _optns = {
            'bkgmod_yield':6,
            'bkghist_yield':hBkg_yield,
            }
        fit_options.update(_optns)

    ptCutTag = args.ptCut
    ptCutProbe = args.ptCut

    outDir = args.dirOut if args.dirOut.endswith("/") else args.dirOut+"/"
    if not os.path.isdir(outDir):
        os.mkdir(outDir)

    outCSVDir = outDir+"csvFiles/"
    if not os.path.isdir(outCSVDir):
        try:
            os.mkdir(outCSVDir)
        except OSError:
            log.warning("directory already exists ...")

    ########################################

    with open(args.mcCorrections) as json_file:
        corr = json.load(json_file)

    ########### Constant settings ##########
    secPerLS = float(23.3)
    currentYear = 2017

    #from 2017H low PU, w/o PU corrections in pb. |eta| < 2.4,  pt>30
    #sigma_fid = 604.280260378
    sigma_fid = 608.305132973

    MassBin_ = 60
    MassMin_ = 56.
    MassMax_ = 116.

    maximumLS = 5000
    LSperMeasurement = 100  # required number of lumi sections per measurement
    LumiPerMeasurement = 20  # minimum recorded lumi for one measurement in pb-1

    ZBBAcc = 0.077904  # Fraction of Z events where both muons are in barrel region
    ZBEAcc = 0.117200  # barrel and endcap
    ZEEAcc = 0.105541  # both endcap

    fullAcc = (ZBBAcc + ZBEAcc + ZEEAcc)

    ZBBAcc = ZBBAcc / fullAcc
    ZBEAcc = ZBEAcc / fullAcc
    ZEEAcc = ZEEAcc / fullAcc

    #nonfiguration of fit models

    log.info("Loading C marco...")
    ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(
        __file__)) + "/calculateDataEfficiency.C")  # load function getZyield(...) and calculateDataEfficiency(...)


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
        hYield = ROOT.TH1D("h_mass_yield_Z","",MassBin_, MassMin_, MassMax_)
        hHLTpassB = ROOT.TH1D("h_mass_HLT_pass_central","",MassBin_, MassMin_, MassMax_)
        hHLTfailB = ROOT.TH1D("h_mass_HLT_fail_central","",MassBin_, MassMin_, MassMax_)
        hHLTpassE = ROOT.TH1D("h_mass_HLT_pass_forward","",MassBin_, MassMin_, MassMax_)
        hHLTfailE = ROOT.TH1D("h_mass_HLT_fail_forward","",MassBin_, MassMin_, MassMax_)
        hSITpassB = ROOT.TH1D("h_mass_SIT_pass_central","",MassBin_, MassMin_, MassMax_)
        hSITfailB = ROOT.TH1D("h_mass_SIT_fail_central","",MassBin_, MassMin_, MassMax_)
        hSITpassE = ROOT.TH1D("h_mass_SIT_pass_forward","",MassBin_, MassMin_, MassMax_)
        hSITfailE = ROOT.TH1D("h_mass_SIT_fail_forward","",MassBin_, MassMin_, MassMax_)
        hGlopassB = ROOT.TH1D("h_mass_Glo_pass_central","",MassBin_, MassMin_, MassMax_)
        hGlofailB = ROOT.TH1D("h_mass_Glo_fail_central","",MassBin_, MassMin_, MassMax_)
        hGlopassE = ROOT.TH1D("h_mass_Glo_pass_forward","",MassBin_, MassMin_, MassMax_)
        hGlofailE = ROOT.TH1D("h_mass_Glo_fail_forward","",MassBin_, MassMin_, MassMax_)

        for run, data_run in byLS_data.groupby('run'):
            if run < int(args.beginRun) or run >= int(args.endRun):
                continue

            LSlist = data_run.query('ls <= {0}'.format(maximumLS))['ls'].values.tolist()

            # check if run was processed already
            log.info("===Loading input DQMIO.root file...")
            eosFileList = glob.glob(eosDir + '/*/*' + str(run) + '*root')
            if not len(eosFileList) > 0:
                log.info("The file does not yet exist for run: " + str(run))
                continue
            else:
                eosFile = eosFileList[0]

            recLumi += sum(data_run.loc[data_run['ls'].isin(LSlist)]['recorded(/pb)'].values)
            dqmfile = ROOT.TFile(eosFile)
            def load(name_):
                return load_histo(name_, dqmfile, LSlist, prefix_dqm, run)

            hYield.Add(load("h_mass_yield_Z"))
            hHLTpassB.Add(load("h_mass_HLT_pass_central"))
            hHLTfailB.Add(load("h_mass_HLT_fail_central"))
            hHLTpassE.Add(load("h_mass_HLT_pass_forward"))
            hHLTfailE.Add(load("h_mass_HLT_fail_forward"))
            hSITpassB.Add(load("h_mass_SIT_pass_central"))
            hSITfailB.Add(load("h_mass_SIT_fail_central"))
            hSITpassE.Add(load("h_mass_SIT_pass_forward"))
            hSITfailE.Add(load("h_mass_SIT_fail_forward"))
            hGlopassB.Add(load("h_mass_Glo_pass_central"))
            hGlofailB.Add(load("h_mass_Glo_fail_central"))
            hGlopassE.Add(load("h_mass_Glo_pass_forward"))
            hGlofailE.Add(load("h_mass_Glo_fail_forward"))

        outSubDir = outDir + "Fits_Inclusive/"
        os.mkdir(outSubDir)

        result = measurement(recLumi, 0, outSubDir, hYield,
            hHLTpassB, hHLTfailB, hHLTpassE, hHLTfailE,
            hSITpassB, hSITfailB, hSITpassE, hSITfailE,
            hGlopassB, hGlofailB, hGlopassE, hGlofailE,
            **fit_options
            )

    log.info("===Looping over runs...")
    for run, data_run in byLS_data.groupby('run'):
        if run < int(args.beginRun) or run >= int(args.endRun):
            continue

        fill = data_run.drop_duplicates('fill')['fill'].values[0]
        LSlist = data_run.query('ls <= {0}'.format(maximumLS))['ls'].values.tolist()

        # check if run was processed already
        outSubDir = outDir + "Run{0}/".format(run)
        if os.path.isdir(outSubDir):
            log.info("Run %i was already processed, skipping and going to next run", run)
            continue
        os.mkdir(outSubDir)

        log.info("===Running Run %i", run)
        log.info("===Running Fill %i", fill)


        log.info("===Loading input DQMIO.root file...")
        eosFileList = glob.glob(eosDir + '/*/*' + str(run) + '*root')

        if not len(eosFileList) > 0:
            log.info("The file does not yet exist for run: " + str(run))
            continue
        else:
            eosFile = eosFileList[0]


        log.info("The file exists: " + str(eosFile) + " for run  " + str(run))
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
            dqmfile = ROOT.TFile(eosFile)

            def load(name_):
                return load_histo(name_, dqmfile, goodLSlist, prefix_dqm, run)

            if not args.inclusive:
                result = measurement(recLumi, len(results), outSubDir,
                    load("h_mass_yield_Z"),
                    load("h_mass_HLT_pass_central"),
                    load("h_mass_HLT_fail_central"),
                    load("h_mass_HLT_pass_forward"),
                    load("h_mass_HLT_fail_forward"),
                    load("h_mass_SIT_pass_central"),
                    load("h_mass_SIT_fail_central"),
                    load("h_mass_SIT_pass_forward"),
                    load("h_mass_SIT_fail_forward"),
                    load("h_mass_Glo_pass_central"),
                    load("h_mass_Glo_fail_central"),
                    load("h_mass_Glo_pass_forward"),
                    load("h_mass_Glo_fail_forward"),
                    **fit_options
                    )

            df = data_run.loc[data_run['ls'].isin(goodLSlist)]

            hist_Zyield = dqmfile.Get("{0}h_mass_yield_Z".format(prefix_dqm))
            hist_ZyieldBB = dqmfile.Get("{0}h_yieldBB_Z".format(prefix_dqm))
            hist_ZyieldEE = dqmfile.Get("{0}h_yieldEE_Z".format(prefix_dqm))

            result.update(ls_corrections(df, result, len(results), hist_Zyield, hist_ZyieldBB, hist_ZyieldEE))

            results.append(result)

        ## Write per measurement csv file - one per run
        log.info("===Writing per Run CSV file")
        results = pd.concat([pd.DataFrame([result]) for result in results])

        with open(outCSVDir + '/csvfile{0}.csv'.format(run), 'w') as file:
            results.to_csv(file, index=False)


    if args.writeSummaryCSV:
        writeSummaryCSV(outCSVDir)

    log.info("===Done")
