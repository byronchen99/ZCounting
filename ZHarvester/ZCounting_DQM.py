#!/usr/bin/env python3

import ROOT
import pandas as pd
import glob
import os
import pdb
import datetime
import numpy as np

from python import utils, logging, common

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None

# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

# -------------------------------------------------------------------------------------
def extract_results(directory, m, cHLT):
    log.info("Extracting fit results in {0} for {1}".format(directory,m))  

    # --- For identification (ID) efficiency        
    NsigHLT2, chi2HLT2 = utils.open_workspace_yield(directory, "HLT_{0}_2".format(etaRegion), m)
    NsigHLT1, chi2HLT1 = utils.open_workspace_yield(directory, "HLT_{0}_1".format(etaRegion), m)
    NsigIDFail, chi2ID = utils.open_workspace_yield(directory, "Sel_{0}_0".format(etaRegion), m)

    NsigHLT2 = utils.unorm(NsigHLT2)
    NsigHLT1 = utils.unorm(NsigHLT1)
    NsigIDFail = utils.unorm(NsigIDFail)

    effHLT = (2 * NsigHLT2) / (2 * NsigHLT2 + NsigHLT1)
    effID = (2 * NsigHLT2 + NsigHLT1) / (2 * NsigHLT2 + NsigHLT1 + NsigIDFail)

    NsigGloFail, chi2GloFail = utils.open_workspace_yield(directory, "Glo_{0}_0".format(etaRegion), m)

    NsigGloFail = utils.unorm(NsigGloFail)

    effGlo = (2 * NsigHLT2 + NsigHLT1 + NsigIDFail) / (2 * NsigHLT2 + NsigHLT1 + NsigIDFail + NsigGloFail)

    effMu = effID*effGlo

    res = {
        "NsigHLT2": NsigHLT2,
        "NsigHLT1": NsigHLT1,
        "NsigIDFail": NsigIDFail,
        "NsigGloFail": NsigGloFail,
        "chi2HLT2": chi2HLT2,
        "chi2HLT1": chi2HLT1,
        "chi2GloFail": chi2GloFail,
        "effHLT": effHLT,
        "effID": effID,
        "effGlo": effGlo,
        "effMu": effMu,
        "cHLT": cHLT,
        "recZCount": (NsigHLT2 + 0.5*NsigHLT1)**2/NsigHLT2 / effMu**2 * cHLT
    }

    return res

################################################################################
if __name__ == '__main__':

    parser = common.parser()
    parser = common.parser_zharvest(parser)

    parser.add_argument("-m", "--measurement", type=int, default=None, 
        help="Only fit a specific measurement of the run")

    common.set_parser_default(parser, "mass", [66,116,50])
    common.set_parser_default(parser, "ptCut", 27)     

    args = parser.parse_args()

    log = logging.setup_logger(__file__, args.verbose)


    ########################################
    # link to resouces
    prefix_dqm="ZCountingInOut-V17_38-"
    cmsswbase = os.environ['CMSSW_BASE']
    resPath = cmsswbase + "/src/ZCounting/ZHarvester/res/"
    if( args.beginRun >= 272007 and args.beginRun < 278808
        # there is an overlap for 2016 F in runs with pre and post VFP settings
        and args.beginRun not in [278769, 278801, 278802, 278803, 278804, 278805, 278808]
    ):                                                          # 2016 pre VFP
        year = 2016
        byLsCSV          = resPath+"/FillByLs_2016.csv"
        mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2016preVFP.root"
        sigTemplates     = args.input+"/"+prefix_dqm+"Summer16preVFP-DYJetsToLL_M_50_LO.root"
    elif args.beginRun < 294645:                                # 2016 post VFP
        year = 2016
        byLsCSV          = resPath+"/FillByLs_2016.csv"
        mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2016postVFP.root"
        sigTemplates     = args.input+"/"+prefix_dqm+"Summer16postVFP-DYJetsToLL_M_50_LO.root"
    elif args.beginRun > 297020 and args.beginRun < 306828:     # 2017
        year = 2017
        byLsCSV          = resPath+"/FillByLs_2017_IsoMu24.csv"
        mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2017.root"
        sigTemplates     = args.input+"/"+prefix_dqm+"Fall17-DYJetsToLL_M_50_LO.root"
    elif args.beginRun >= 306926 and args.beginRun < 307083:    # 2017 H
        year = 2017
        byLsCSV          = resPath+"/FillByLs_2017_lowPU.csv"
        mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2017.root"
        sigTemplates     = args.input+"/"+prefix_dqm+"Fall17-DYJetsToLL_M_50_LO.root"
    elif args.beginRun >= 315252 and args.beginRun < 325273:    # 2018
        year = 2018
        byLsCSV          = resPath+"/FillByLs_2018.csv"
        mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2018.root"
        sigTemplates     = args.input+"/"+prefix_dqm+"Autumn18-DYJetsToLL_M_50_LO.root"
    elif args.beginRun >= 355100 and args.beginRun < 362760:    # 2022
        year = 2022
        byLsCSV = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/brilcalcByLS/byLS_Collisions22_355100_362760_Muon_20230210.csv"
        mcCorrelations   = f"{args.input}/2022/CorrelationFactors/c_nPV_2022.root"
        prefix_dqm =  "ZCountingAll-V01-"
        sigTemplates = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/SignalTemplates/ZCountingAll-V01-Winter22-DYJetsToLL_M_50_LO.root"
    elif args.beginRun >= 362760:                               # 2023
        year = 2023
    else:
        log.warning("specified begin run {0} unknown! exit()".format(args.beginRun))
        exit()

    energy = 13.6 if year >= 2022 else 13

    byLsCSV = byLsCSV          if args.byLsCSV       == "default"   else args.byLsCSV
    eosDir  = args.input

    measurement      = args.measurement

    log.info("----------------------------------")
    log.info("Use eosDir:              {0}".format(eosDir))
    log.info("Use byLsCSV:             {0}".format(byLsCSV))
    log.info("Use sigTemplates:        {0}".format(sigTemplates))
    log.info("Mass range from:         {0} to {1}".format(*args.mass))
    log.info("Lumi per Measurement:    {0}".format(args.LumiPerMeasurement))
    log.info("----------------------------------")
    
    sigModel = common.sigModels[args.sigModel]
    sigModelSta = common.sigModels[args.sigModel.replace("Gen","MCxGauss")]    # for standalone we don't use Gen
    bkgModelPass = common.bkgModelsPass[args.bkgModel]
    bkgModelFail = common.bkgModelsFail[args.bkgModel]

    ########## Input configuration ##########
    # ByLS csv inputs generated by testBRIL.sh
    byLS_filelist = glob.glob(byLsCSV)
    byLS_filelist.sort(key=os.path.getmtime)
    byLS_filename = byLS_filelist[-1]
    log.info("The brilcalc csv file: " + str(byLS_filename))

    outDir = args.output if args.output.endswith("/") else args.output+"/"
    if not os.path.isdir(outDir):
        os.mkdir(outDir)

    outCSVDir = outDir+"csvFiles/"
    if not os.path.isdir(outCSVDir):
        try:
            os.mkdir(outCSVDir)
        except OSError:
            log.warning("Directory already exists ...")

    ########### Constant settings ##########
    secPerLS = float(23.3)
    
    LumiPerMeasurement = args.LumiPerMeasurement  # minimum recorded lumi for one measurement in pb-1

    #configuration of fit models
    MassMin_ = int(args.mass[0])
    MassMax_ = int(args.mass[1])
    MassBin_ = int(args.mass[2])
    MassBinWidth = (MassMax_ - MassMin_)/MassBin_

    MassMinSta_ = int(50)
    MassMaxSta_ = int(130)
    MassBinSta_ = int((MassMaxSta_ - MassMinSta_)/MassBinWidth)

    npvMin_ = 0.5     
    npvMax_ = 100.5   

    npvBin_ = int(npvMax_ - npvMin_)

    maximumLS = 2500        # maximum luminosity block number that is stored in the DQM histograms

    if args.etaCut == 0.9:
        etaRegion = "B"
        etaRegionZ = "BB"
    else:
        etaRegion = "I"
        etaRegionZ = "I"

    if not args.collect:
        log.info("Loading C marcos...")
        # load functions for fitting
        ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/Utils/calculateDataEfficiency.C")

        if args.bkgModel == "Das":
            ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/Utils/RooGaussDoubleSidedExp.cc+")

        ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

        ROOT.set_npvRange(npvMin_, npvMax_)
        if year >= 2022:
            ROOT.set_energy(energy)

        ROOT.set_ptCut(args.ptCut)
        ROOT.set_etaCut(args.etaCut)
    
    byLS_data = utils.load_input_csv(byLS_filename)
    byLS_data = byLS_data.loc[(byLS_data['run'] >= int(args.beginRun)) & (byLS_data['run'] < int(args.endRun))]
    
    #####################################   
    recLumi = 0
    results = []
    log.info(f"Looping over runs... {args.beginRun} to {args.endRun}")
    for run, byLS_run in byLS_data.groupby('run', sort=True):

        fill = byLS_run.drop_duplicates('fill')['fill'].values[0]
        LSlist = byLS_run['ls'].values.tolist()
        Lumilist = byLS_run.loc[byLS_run['ls'].isin(LSlist)]['recorded(/pb)'].values.tolist()

        log.info(f"Running Fill {fill}")
        log.info(f"Running Run {run}")

        log.info(f"Now at run {run}")
        fileName = utils.getFileName(eosDir,run)
        if fileName is None:
            continue
        log.info(f"Found file `{fileName}`")

        log.debug("Have lumi secion list {0}".format(LSlist))        
        log.info("Looping over measurements...")
        for m, goodLSlist in enumerate(
            utils.get_ls_for_next_measurement(lumisections=LSlist, luminosities=Lumilist, #zcounts=ZCountlist, 
                lumiPerMeasurement=LumiPerMeasurement)
        ):
            if measurement is not None and measurement < m:
                break
                    
            # create datafram byLS for measurement
            byLS_m = byLS_run.loc[byLS_run['ls'].isin(goodLSlist)]

            # get histograms from good lumisections
            log.info("Load histograms ...")
            prefix=f"DQMData/Run {run}/ZCounting/Run summary/Histograms/"
            
            # get histogram with primary vertex distribution
            hPV = utils.load_histogram("h_npv", fileName, goodLSlist, run=run, prefix=prefix, suffix="new", pileup=True)

            # get histograms binned in mass
            def load(name_):
                return utils.load_histogram(name_, fileName, goodLSlist, run=run, 
                    MassBin=MassBin_, MassMin=MassMin_, MassMax=MassMax_, 
                    prefix=prefix, 
                    suffix="new")
            
            # load histograms for hlt efficiency and Z yield
            h2HLT = load(f"h_mass_2HLT_BB")
            h1HLT = load(f"h_mass_1HLT_BB")

            # load histograms with probes that fail selection, for selection efficiency
            hIDfail = load(f"h_mass_SIT_fail_BB")

            # load histograms for global muon efficiency
            # hGlopass = load(f"h_mass_Glo_pass_BB")
            hGlofail = load(f"h_mass_Glo_fail_BB")

            # load histograms for standalone muon efficiency
            # hStapass = load(f"h_mass_Sta_pass_BB")
            # hStafail = load(f"h_mass_Sta_fail_BB")

            if etaRegion == "I":
                h2HLT += load(f"h_mass_2HLT_BE")
                h2HLT += load(f"h_mass_2HLT_EE")

                h1HLT += load(f"h_mass_1HLT_BE")
                h1HLT += load(f"h_mass_1HLT_EE")

                hIDfail += load(f"h_mass_SIT_fail_BE")
                hIDfail += load(f"h_mass_SIT_fail_EE")

                # hGlopass += load(f"h_mass_Glo_pass_BE")
                # hGlopass += load(f"h_mass_Glo_pass_EE")
                hGlofail += load(f"h_mass_Glo_fail_BE")
                hGlofail += load(f"h_mass_Glo_fail_EE")
                
                # hStapass += load(f"h_mass_Sta_pass_BE")
                # hStapass += load(f"h_mass_Sta_pass_EE")
                # hStafail += load(f"h_mass_Sta_fail_BE")
                # hStafail += load(f"h_mass_Sta_fail_EE")

            recLumi = byLS_m['recorded(/pb)'].sum()

            log.info("Have now recorded lumi = {0}".format(recLumi))            

            log.info("Histograms filled ...")  

            outSubDir = outDir + "Run{0}/".format(run)

            log.debug("Running measurement {0}".format(m))

            if not args.collect:
                
                if measurement is None or measurement == m:
                    # skip the fit if we look for another measurement
                
                    if not os.path.isdir(outSubDir):
                        os.mkdir(outSubDir)
                    
                    ROOT.set_output(outSubDir)
                    ROOT.set_luminosity(recLumi)

                    hPV = utils.np_to_hist(hPV, npvBin_, npvMin_, npvMax_, "TH1", "hPV")

                    def to_hist(h_, name):
                        return utils.np_to_hist(h_, MassBin_, MassMin_, MassMax_, "TH1", name)
                    
                    h2HLT = to_hist(h2HLT, "HLT2")
                    h1HLT = to_hist(h1HLT, "HLT1")
                    hIDfail = to_hist(hIDfail, "IDfail")
                    # hGlopass = to_hist(hGlopass, "Glopass")
                    hGlofail = to_hist(hGlofail, "Glofail")
                    # hStapass = to_hist(hStapass, "Stapass")
                    # hStafail = to_hist(hStafail, "Stafail")

                    ROOT.getZyield(h2HLT, m, "HLT", etaRegion, sigModel, bkgModelPass, 2, sigTemplates, 0)
                    ROOT.getZyield(h1HLT, m, "HLT", etaRegion, sigModel, bkgModelPass, 1, sigTemplates, 0)
                    ROOT.getZyield(hIDfail, m, "Sel", etaRegion, sigModel, bkgModelPass, 0, sigTemplates, 0)

                    ROOT.set_massRange(MassMinSta_, MassMaxSta_, MassBinSta_)
                    # ROOT.getZyield(hGlopass, m, "Glo", etaRegion, sigModel, bkgModelPass, 1, sigTemplates, 0)
                    ROOT.getZyield(hGlofail, m, "Glo", etaRegion, sigModel, bkgModelFail, 0, sigTemplates, 0)
                    ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

                    # ROOT.getZyield(hStapass, m, "Sta", etaRegion, sigModel, bkgModelPass, 1, sigTemplates, 0)
                    # ROOT.getZyield(hStafail, m, "Sta", etaRegion, sigModel, bkgModelFail, 0, sigTemplates, 0)

                    # remove the histogram templates, not needed anymore
                    os.system("rm {0}/histTemplates_*".format(outSubDir))

            cHLT = ROOT.extractCorrelation_HLT(sigTemplates, hPV, etaRegionZ)

            result = extract_results(outSubDir, m, cHLT)

            if result:
                delLumi = byLS_m['delivered(/pb)'].sum()
                recLumi = byLS_m['recorded(/pb)'].sum()

                # deadtime - fraction during the good lumisections
                deadtime = recLumi / delLumi
                timewindow = len(byLS_m) * secPerLS

                # convert time string to datetime format
                beginTime = utils.to_DateTime(byLS_m['time'][0], string_format = "mm/dd/yy")
                endTime = utils.to_DateTime(byLS_m['time'][-1], string_format = "mm/dd/yy") + datetime.timedelta(seconds=secPerLS)
                
                # total time window from the beginning of the first to the end of the last lumisection
                totaltimewindow = (endTime - beginTime).total_seconds()

                # delivered, efficiency and deadtime corrected, Z boson count
                delZCount = result["recZCount"] / deadtime * (totaltimewindow / timewindow)

                result.update({
                    "fill": fill,
                    "run": run,
                    "measurement": m,
                    "beginTime": beginTime.strftime("%y/%m/%d %H:%M:%S"),
                    "endTime": endTime.strftime("%y/%m/%d %H:%M:%S"),
                    "delLumi": delLumi,
                    "recLumi": recLumi,
                    "timewindow": timewindow,
                    "pileUp": byLS_m['avgpu'].mean(),
                    "delZCount": delZCount,
                })
            
                results.append(result)
            else:
                log.info("No result - continue")
        
        if measurement is None or measurement == m:
            utils.writePerRunCSV(results, outCSVDir, run)

        results = []

    if args.writeSummaryCSV:
        utils.writeSummaryCSV(outCSVDir)

    log.info("Done")
