
def measurement(recLumi, m, outputDir,
    h1RecoBB, h1RecoBE, h1RecoEE,
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
    # compute Z yield
    ZyieldresBB_m = ROOT.getZyield(h1RecoBB, outputDir, m, sigmod_yield, bkgmod_yield,
        ptCutTag, ptCutProbe, "BB", recLumi, sigtempl_yield, bkghist_yield)

    ZyieldresBE_m = ROOT.getZyield(h1RecoBE, outputDir, m, sigmod_yield, bkgmod_yield,
        ptCutTag, ptCutProbe, "BE", recLumi, sigtempl_yield, bkghist_yield)

    ZyieldresEE_m = ROOT.getZyield(h1RecoEE, outputDir, m, sigmod_yield, bkgmod_yield,
        ptCutTag, ptCutProbe, "EE", recLumi, sigtempl_yield, bkghist_yield)

    # check if there is was enough statistics
    if np.isnan(ZyieldresBB_m[0]) or np.isnan(ZyieldresBE_m[0]) or np.isnan(ZyieldresEE_m[0]) or ZyieldresBB_m[0] * ZyieldresBE_m[0] * ZyieldresEE_m[0] == 0:
        return None

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
    ZBBeff = (GloeffB_m * GloeffB_m * SeleffB_m * SeleffB_m * (1 - (1 - HLTeffB_m) * (1 - HLTeffB_m)))
    ZBEeff = (GloeffB_m * GloeffE_m * SeleffB_m * SeleffE_m * (1 - (1 - HLTeffB_m) * (1 - HLTeffE_m)))
    ZEEeff = (GloeffE_m * GloeffE_m * SeleffE_m * SeleffE_m * (1 - (1 - HLTeffE_m) * (1 - HLTeffE_m)))

    # Statistic Uncertainties (low,high) error propagation
    ZBBeff_EStat = [0., 0.]
    ZBEeff_EStat = [0., 0.]
    ZEEeff_EStat = [0., 0.]
    for i in (1, 2):
        ZBBeff_EStat[i - 1] = 2 * ZBBeff * np.sqrt(
            (GloeffresB_m[i] / GloeffB_m) ** 2 +
            (SeleffresB_m[i] / SeleffB_m) ** 2 +
            ((1 - HLTeffB_m) / (1 - (1 - HLTeffB_m) ** 2) * HLTeffresB_m[i]) ** 2
        )
        ZEEeff_EStat[i - 1] = 2 * ZEEeff * np.sqrt(
            (GloeffresE_m[i] / GloeffE_m) ** 2 +
            (SeleffresE_m[i] / SeleffE_m) ** 2 +
            ((1 - HLTeffE_m) / (1 - (1 - HLTeffE_m) ** 2) * HLTeffresE_m[i]) ** 2
        )
        ZBEeff_EStat[i - 1] = ZBEeff * np.sqrt(
            (GloeffresB_m[i] / GloeffB_m) ** 2 +
            (GloeffresE_m[i] / GloeffE_m) ** 2 +
            (SeleffresB_m[i] / SeleffB_m) ** 2 +
            (SeleffresE_m[i] / SeleffE_m) ** 2 +
            ((1 - HLTeffE_m) / (1 - (1 - HLTeffB_m) * (1 - HLTeffE_m)) * HLTeffresB_m[i]) ** 2 +
            ((1 - HLTeffB_m) / (1 - (1 - HLTeffB_m) * (1 - HLTeffE_m)) * HLTeffresE_m[i]) ** 2
        )


    res = {
        "zYieldBB": ZyieldresBB_m[0],
        "zYieldBB_err": max(ZyieldresBB_m[1:3]),
        "zYieldBB_chi2": ZyieldresBB_m[3],
        "zYieldBB_purity": ZyieldresBB_m[4],
        "zYieldBB_purity_err": max(ZyieldresBB_m[5:]),
        "zYieldBE": ZyieldresBE_m[0],
        "zYieldBE_err": max(ZyieldresBE_m[1:3]),
        "zYieldBE_chi2": ZyieldresBE_m[3],
        "zYieldBE_purity": ZyieldresBE_m[4],
        "zYieldBE_purity_err": max(ZyieldresBE_m[5:]),
        "zYieldEE": ZyieldresEE_m[0],
        "zYieldEE_err": max(ZyieldresEE_m[1:3]),
        "zYieldEE_chi2": ZyieldresEE_m[3],
        "zYieldEE_purity": ZyieldresEE_m[4],
        "zYieldEE_purity_err": max(ZyieldresEE_m[5:]),
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
        "ZBBeff": ZBBeff,
        "ZBEeff": ZBEeff,
        "ZEEeff": ZEEeff,
        "ZBBeff_stat": max(ZBBeff_EStat),
        "ZBEeff_stat": max(ZBEeff_EStat),
        "ZEEeff_stat": max(ZEEeff_EStat),
    }
    res["zBB_relstat"] = res['zYieldBB_err']/ZyieldresBB_m[0] + res["ZBBeff_stat"]/res["ZBBeff"]
    res["zBE_relstat"] = res['zYieldBE_err']/ZyieldresBE_m[0] + res["ZBEeff_stat"]/res["ZBEeff"]
    res["zEE_relstat"] = res['zYieldEE_err']/ZyieldresEE_m[0] + res["ZEEeff_stat"]/res["ZEEeff"]
    return res

def ls_corrections(data_m, result_m, m, h2ZyieldBB, h2ZyieldBE, h2ZyieldEE):
    # --- per ls dataframe, one per measurement

    for eta, hist in (("BB",h2ZyieldBB), ("BE",h2ZyieldBE), ("EE",h2ZyieldEE)):

        data_m["eff{0}_mc".format(eta)] = result_m["Z{0}eff".format(eta)] - (data_m['avgpu'] * corr['{0}_a'.format(eta)] + corr['{0}_b'.format(eta)])

        data_m['yield{0}'.format(eta)] = data_m['ls'].apply(lambda ls: hist.ProjectionY("", ls, ls, "e").GetEntries())

        data_m['zYield{0}'.format(eta)]  = data_m['yield{0}'.format(eta)] * result_m["zYield{0}_purity".format(eta)]
        data_m['zDel{0}'.format(eta)]    = data_m['zYield{0}'.format(eta)] / result_m["Z{0}eff".format(eta)]
        data_m['zDel{0}_mc'.format(eta)] = data_m['zYield{0}'.format(eta)] / data_m['eff{0}_mc'.format(eta)]

        data_m['z{0}_relstat'.format(eta)] = (1./np.sqrt(data_m['zYield{0}'.format(eta)])).replace(np.inf, 0.) \
            + result_m["Z{0}eff_stat".format(eta)]/result_m["Z{0}eff".format(eta)] \
            + result_m["zYield{0}_purity_err".format(eta)]/result_m["zYield{0}_purity".format(eta)]

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
        "ZBBeff_mc": data_m['effBB_mc'].mean(),
        "ZBEeff_mc": data_m['effBE_mc'].mean(),
        "ZEEeff_mc": data_m['effEE_mc'].mean(),
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

    os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
    from ZUtils.python.utils import to_RootTime

    # disable panda warnings when assigning a new column in the dataframe
    pd.options.mode.chained_assignment = None
    # turn off graphical output on screen
    ROOT.gROOT.SetBatch(True)
    cmsswbase = os.environ['CMSSW_BASE']

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--beginRun", help="first run to analyze [%(default)s]", default=272007)
    parser.add_argument("-e", "--endRun", help="analyze stops when comes to this run [%(default)s]", default=1000000)
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

    MassBin_ = 60
    MassMin_ = 56.
    MassMax_ = 116.

    maximumLS = 5000
    LSperMeasurement = 100  # required number of lumi sections per measurement
    LumiPerMeasurement = 20  # minimum recorded lumi for one measurement in pb-1

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
        hYieldBB = ROOT.TH1D("h_mass_yielBBd_Z","",MassBin_, MassMin_, MassMax_)
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
            hGlopassB.Add(load("h_mass_Glo_pass_central"))
            hGlofailB.Add(load("h_mass_Glo_fail_central"))
            hGlopassE.Add(load("h_mass_Glo_pass_forward"))
            hGlofailE.Add(load("h_mass_Glo_fail_forward"))

        outSubDir = outDir + "Fits_Inclusive/"
        os.mkdir(outSubDir)

        result = measurement(recLumi, 0, outSubDir,
            hYieldBB, hYieldBE, hYieldEE,
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

        if not args.inclusive:
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
                    load("h_mass_yieldBB_Z"),
                    load("h_mass_yieldBE_Z"),
                    load("h_mass_yieldEE_Z"),
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

            if result is None:
                continue

            df = data_run.loc[data_run['ls'].isin(goodLSlist)]

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
