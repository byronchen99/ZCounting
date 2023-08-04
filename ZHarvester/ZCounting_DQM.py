#!/usr/bin/env python3

#import ROOT
import pandas as pd
import glob
import os
import pdb
import datetime
import numpy as np
import uncertainties as unc
# from roofit.utils import load_makros
import h5py
from zfit_tools.fit import fit, fit_2
from zfit_tools.get_acceptance import get_acceptance, acceptance_uncertainties
# from templates.fitting_templates import templates_2d_to_1d

from common import utils, logging, parsing

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None

# turn off graphical output on screen
# ROOT.gROOT.SetBatch(True)

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------


def extract_results(data_histograms, templates_1d, correlation_templates, acceptance_templates, hPV):

    # MassBin_, MassMin_, MassMax_ = massbins["Other_Categories"]
    # MassBinSta_, MassMinSta_, MassMaxSta_ = massbins["Glo"]
    data_hist_HLT2 = data_histograms["HLT2"]
    data_hist_HLT1 = data_histograms["HLT1"]
    data_hist_IDfail = data_histograms["IDfail"]
    data_hist_Glopass = data_histograms["Glopass"]
    data_hist_Glofail = data_histograms["Glofail"]
    data_hist_Stapass = data_histograms["Stapass"]
    data_hist_Stafail = data_histograms["Stafail"]
    

    # fitting
    NsigHLT2, alphaHLT2, betaHLT2, NbkgHLT2 = fit(templates_1d, data_hist_HLT2, "HLT2")
    NsigHLT1, alphaHLT1, betaHLT1, NbkgHLT1  = fit(templates_1d, data_hist_HLT1, "HLT1")
    NsigIDFail, alphaIDFail, betaIDFail, NbkgIDFail  = fit(templates_1d, data_hist_IDfail, "IDfail")

    if NsigHLT2 + NsigHLT1 <= 0:
        effHLT = utils.unorm(0)
        effID = utils.unorm(0)
    else:
        effHLT = (2 * NsigHLT2) / (2 * NsigHLT2 + NsigHLT1)
        effID = (2 * NsigHLT2 + NsigHLT1) / (2 * NsigHLT2 + NsigHLT1 + NsigIDFail)


    NsigGloPass, NbkgGloPass = fit_2(templates_1d, data_hist_Glopass, "Glopass") # for now use only a template yield fit as there are no corrections for Glo yet
    NsigGloFail, NbkgGloFail = fit_2(templates_1d, data_hist_Glofail, "Glofail")
    # NsigGloPass, alphaGloPass, betaGloPass, NbkgGloPass = fit(templates_1d, data_hist_Glopass, "Glopass")
    # NsigGloFail, alphaGloFail, betaGloFail, NbkgGloFail  = fit(templates_1d, data_hist_Glofail, "Glofail")
    NsigStaPass, alphaStaPass, betaStaPass, NbkgStaPass = fit(templates_1d, data_hist_Stapass, "Stapass")
    NsigStaFail, alphaStaFail, betaStaFail, NbkgStaFail  = fit(templates_1d, data_hist_Stafail, "Stafail")


    if NsigGloPass + NsigGloFail <= 0:
        effGlo = utils.unorm(0)
    else:
        effGlo = NsigGloPass / (NsigGloPass + NsigGloFail)

    if NsigStaPass + NsigStaFail <= 0:
        effSta = utils.unorm(0)
    else:
        effSta = NsigStaPass / (NsigStaPass + NsigStaFail)

    effMu = effID*effGlo*effSta
    
    
    nPV = hPV.view()
    
    cHLT = sum(correlation_templates["HLT"] * nPV) / sum(nPV)
    cID = sum(correlation_templates["ID_total"] * nPV) / sum(nPV)
    cIO = sum(correlation_templates["io"] * nPV) / sum(nPV)

    
    # get acceptance
    morphing_params = morphing_params = {"HLT2_alpha": alphaHLT2, "HLT1_alpha": alphaHLT1, "IDfail_alpha": alphaIDFail, "Stapass_alpha": alphaStaPass, "Stafail_alpha": alphaStaFail,
     "HLT2_beta": betaHLT2, "HLT1_beta": betaHLT1, "IDfail_beta": betaIDFail, "Stapass_beta": betaStaPass, "Stafail_beta": betaStaFail}
    
    morphing_params_nominal, morphing_params_unc = {}, {}
    for key in morphing_params.keys():
        morphing_params_nominal[key] = morphing_params[key].nominal_value
        morphing_params_unc[key] = morphing_params[key].std_dev
   
    
    acceptance_by_PV = get_acceptance(morphing_params_nominal, acceptance_templates, plot = True)
    acceptance_0 = sum(acceptance_by_PV * nPV) / sum(nPV)

    acceptance_uncertainty_plus_PV, acceptance_uncertainty_minus_PV = acceptance_uncertainties(morphing_params_nominal, morphing_params_unc, acceptance_templates, acceptance_by_PV)
    acceptance_uncertainty_plus = sum(acceptance_uncertainty_plus_PV * nPV) / sum(nPV)
    acceptance_uncertainty_minus = sum(acceptance_uncertainty_minus_PV * nPV) / sum(nPV)
    acceptance_uncertainty = max(acceptance_uncertainty_plus, acceptance_uncertainty_minus)

    acceptance = unc.ufloat(acceptance_0, acceptance_uncertainty)


    if effMu <= 0 or NsigHLT2 <=0:
        recZCount = utils.unorm(0)
    else:
        recZCount = (NsigHLT2 + 0.5*NsigHLT1)**2/NsigHLT2 / effMu**2 * cHLT * cID / cIO**2 * acceptance  # Note the divide by cIO!
        

    res = {
    # Signal numbers in mutual exclusive categories
    "NsigHLT2": NsigHLT2,
    "NsigHLT1": NsigHLT1,
    "NsigIDFail": NsigIDFail,
    "NsigGloPass": NsigGloPass,
    "NsigGloFail": NsigGloFail,
    "NsigStaPass": NsigStaPass,
    "NsigStaFail": NsigStaFail,
    # Background numbers in mutual exclusive categories
    "NbkgHLT1": NbkgHLT1,
    "NbkgHLT2": NbkgHLT2,
    "NbkgIDFail": NbkgIDFail,
    "NbkgGloPass": NbkgGloPass,
    "NbkgGloFail": NbkgGloFail,
    "NbkgStaPass": NbkgStaPass,
    "NbkgStaFail": NbkgStaFail,
    # Chi2s from fits
    # "chi2HLT2": chi2HLT2,
    # "chi2HLT1": chi2HLT1,
    # "chi2GloPass": chi2GloPass,
    # "chi2GloFail": chi2GloFail,
    # "chi2StaPass": chi2StaPass,
    # "chi2StaFail": chi2StaFail,
    # Efficiencies
    "effHLT": effHLT,
    "effID": effID,
    "effGlo": effGlo,
    "effSta": effSta,
    "effMu": effMu,
    # Corrections from MC
    "cHLT": cHLT,
    "cIO": cIO,
    "cID": cID,
    # "cAcceptanceInner": cAcceptanceInner,
    # "cAcceptanceOuter": cAcceptanceOuter,
    "acceptance": acceptance,
    # Efficiency corrected number of Z bosons
    "recZCount": recZCount,  
    # # Fitting Parameters
    "alphaHLT2": alphaHLT2,
    "betaHLT2": betaHLT2,
    "alphaHLT1": alphaHLT1,
    "betaHLT1": betaHLT1,
    "alphaIDFail": alphaIDFail,
    "betaIDFail": betaIDFail,
    # "alphaGloPass": alphaGloPass,
    # "betaGlopass": betaGloPass,
    # "alphaGloFail": alphaGloFail,
    # "betaGloFail": betaGloFail,
    "alphaStaPass": alphaStaPass,
    "betaStaPass": betaStaPass,
    "alphaStaFail": alphaStaFail,
    "betaStaFail": betaStaFail
    }

    return res


################################################################################
if __name__ == '__main__':

    parser = parsing.parser_zharvest()

    parser.add_argument("-m", "--measurement", type=int, default=None, 
        help="Only fit a specific measurement of the run")
    
    parsing.set_parser_default(parser, "mass", [66,116,50])
    parsing.set_parser_default(parser, "ptCut", 27)     

    args = parser.parse_args()

    log = logging.setup_logger(__file__, args.verbose)


    ########################################
    # link to resouces
    year = utils.get_year_by_run(args.beginRun)
    energy = 13.6 if year >= 2022 else 13

    # cmsswbase = os.environ['CMSSW_BASE']
    # resPath = cmsswbase + "/src/ZCounting/ZHarvester/res/"
    # if( args.beginRun >= 272007 and args.beginRun < 278808
    #     # there is an overlap for 2016 F in runs with pre and post VFP settings
    #     and args.beginRun not in [278769, 278801, 278802, 278803, 278804, 278805, 278808]
    # ):                                                          # 2016 pre VFP
    #     byLsCSV          = resPath+"/FillByLs_2016.csv"
    #     mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2016preVFP.root"
    # elif args.beginRun < 294645:                                # 2016 post VFP
    #     byLsCSV          = resPath+"/FillByLs_2016.csv"
    #     mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2016postVFP.root"
    # elif args.beginRun > 297020 and args.beginRun < 306828:     # 2017
    #     byLsCSV          = resPath+"/FillByLs_2017_IsoMu24.csv"
    #     mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2017.root"
    # elif args.beginRun >= 306926 and args.beginRun < 307083:    # 2017 H
    #     byLsCSV          = resPath+"/FillByLs_2017_lowPU.csv"
    #     mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2017.root"
    # elif args.beginRun >= 315252 and args.beginRun < 325273:    # 2018
    #     byLsCSV          = resPath+"/FillByLs_2018.csv"
    #     mcCorrelations   = resPath+"mcCorrections/V17_38/c_nPV_2018.root"
    # elif args.beginRun >= 355100 and args.beginRun < 362760:    # 2022
    #     byLsCSV = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/brilcalcByLS/byLS_Collisions22_355100_362760_Muon_20230210.csv"
    #     mcCorrelations  = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/CorrelationFactors/MCClosure_V19_07/c_nPV_2022.root"
    #     sigTemplates = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/SignalTemplates/ZCountingInOut-V19_07-Winter22-DYJetsToLL_M_50_LO.root"
    # elif args.beginRun >= 366403:    # 2023
    #     byLsCSV = "/eos/cms/store/group/comm_luminosity/ZCounting/2023/brilcalcByLS/byLS_Collisions23_366403_367111_Muon_20230508.csv"
    #     mcCorrelations  = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/CorrelationFactors/MCClosure_V19_07/c_nPV_2022.root"
    #     sigTemplates = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/SignalTemplates/ZCountingInOut-V19_07-Winter22-DYJetsToLL_M_50_LO.root"


    byLsCSV = args.byLsCSV
    eosDir  = args.input
    sigTemplates = None

    measurement      = args.measurement

    log.info("----------------------------------")
    log.info("Use eosDir:              {0}".format(eosDir))
    log.info("Use byLsCSV:             {0}".format(byLsCSV))
    log.info("Use sigTemplates:        {0}".format(sigTemplates))
    log.info("Mass range from:         {0} to {1}".format(*args.mass))
    log.info("Lumi per Measurement:    {0}".format(args.LumiPerMeasurement))
    log.info("----------------------------------")
    
    sigModel = parsing.sigModels[args.sigModel]
    sigModelSta = parsing.sigModels[args.sigModel.replace("Gen","MCxGauss")]    # for standalone we don't use Gen
    bkgModelPass = parsing.bkgModelsPass[args.bkgModel]
    bkgModelFail = parsing.bkgModelsFail[args.bkgModel]

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
    secPerLS = float(23.31112)
    
    LumiPerMeasurement = args.LumiPerMeasurement  # minimum recorded lumi for one measurement in pb-1

    #configuration of fit models
    MassMin_ = int(args.mass[0])
    MassMax_ = int(args.mass[1])
    MassBin_ = int(args.mass[2])
    MassBinWidth = (MassMax_ - MassMin_)/MassBin_

    MassMinSta_ = int(min(MassMin_, 56))
    MassMaxSta_ = int(max(MassMax_, 126))
    MassBinSta_ = int((MassMaxSta_ - MassMinSta_)/MassBinWidth)

    npvMin_ = 0.5     
    npvMax_ = 100.5   

    npvBin_ = int(npvMax_ - npvMin_)

    maximumLS = 2500        # maximum luminosity block number that is stored in the DQM histograms

    if args.etaCut == 0.9 and args.etaMin == 0.0:
        etaRegion = "B"
        etaRegionZ = "BB"
        etaRegions = ["BB"]
    elif args.etaCut == 2.4 and args.etaMin == 0.9:
        etaRegion = "E"
        etaRegionZ = "EE"
        etaRegions = ["EE"]
    elif args.etaCut == 2.4 and args.etaMin == 0.0:
        etaRegion = "I"
        etaRegionZ = "I"
        etaRegions = ["BB","BE","EE"]
    else:
        log.error(f"Eta cut and {args.etaMin} < |eta| < {args.etaCut} not supported for DQM histograms!")

    #etaRegions = ["BB","BE","EE"] if etaRegion == "I" else ["BB"]

    # if not args.collect:        
    #     load_makros(args.bkgModel, args.ptCut, args.etaMin, args.etaCut, year, MassMin_, MassMax_, MassBin_, npvMin_, npvMax_)


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

        # path in root file where the histograms are stored
        dqmPrefix=f"DQMData/Run {run}/ZCounting/Run summary/Histograms/"

        # get list of Z boson candidates per lumisection
        ZCountlist = utils.load_histogram(
            [f"h_mass_2HLT_{r}" for r in etaRegions]+[f"h_mass_1HLT_{r}" for r in etaRegions], 
            fileName, LSlist, run=run, 
            prefix=dqmPrefix, 
            entries=True)

        log.info(f"Looping over measurements in intervals of {LumiPerMeasurement}/pb")
        for m, goodLSlist in enumerate(
            utils.get_ls_for_next_measurement(lumisections=LSlist, 
                #luminosities=Lumilist, 
                # zcounts=ZCountlist, 
                lumiPerMeasurement=LumiPerMeasurement)
        ):
            if measurement is not None and measurement < m:
                break

            log.debug("Running measurement {0}".format(m))
            log.debug("Good lumi list: {0}".format(goodLSlist))

            # create datafram byLS for measurement
            byLS_m = byLS_run.loc[byLS_run['ls'].isin(goodLSlist)]

            # get histograms from good lumisections
            log.info("Load histograms ...")
            
            # get histogram with primary vertex distribution
            hPV = utils.load_histogram("h_npv", fileName, goodLSlist, run=run, prefix=dqmPrefix, pileup=True)
            hPV = utils.np_to_hist(hPV, npvBin_, npvMin_, npvMax_, "hist", "hPV")

            # get histograms binned in mass
            def load(name_, mBins=MassBin_, mMin=MassMin_, mMax=MassMax_):
                hist_np = utils.load_histogram([f"h_mass_{name_}_{r}" for r in etaRegions], fileName, goodLSlist, run=run, 
                    MassBin=mBins, MassMin=mMin, MassMax=mMax, 
                    prefix=dqmPrefix)
                return utils.np_to_hist(hist_np, mBins, mMin, mMax, "hist", name_) #convert numpy array into ROOT.TH1

            # load histograms for hlt efficiency and Z yield
            h2HLT = load("2HLT")
            h1HLT = load("1HLT")

            # load histograms with probes that fail ID, for ID efficiency
            hIDfail = load("ID_fail")

            # load histograms for global muon efficiency
            hGlopass = load("Glo_pass", MassBinSta_, MassMinSta_, MassMaxSta_)
            hGlofail = load("Glo_fail", MassBinSta_, MassMinSta_, MassMaxSta_)

            # load histograms for standalone muon efficiency
            hStapass = load("Sta_pass")
            hStafail = load("Sta_fail")

            recLumi = byLS_m['recorded(/pb)'].sum()

            log.info("Have now recorded lumi = {0}".format(recLumi))            
            log.info("Have now {0} ({1}) HLT 2 (1) events".format(sum(h2HLT), sum(h1HLT)))

            log.info("Histograms filled ...")  

            outSubDir = outDir + "Run{0}/".format(run)

            log.debug("Running measurement {0}".format(m))

            if not args.collect:
                
                if measurement is None or measurement == m:
                    # skip the fit if we look for another measurement
                
                    if not os.path.isdir(outSubDir):
                        os.mkdir(outSubDir)

                    #result = fit(h2HLT, MassBin_, MassMin_, MassMax_, "HLT_2")
                    #pdb.set_trace()

                    # ROOT.set_output(outSubDir)
                    # ROOT.set_luminosity(recLumi)

                    # ROOT.getZyield(h2HLT, m, "HLT", etaRegionZ, sigModel, bkgModelPass, 2, sigTemplates, 0)
                    # ROOT.getZyield(h1HLT, m, "HLT", etaRegionZ, sigModel, bkgModelPass, 1, sigTemplates, 0)
                    # ROOT.getZyield(hIDfail, m, "ID", etaRegion, sigModel, bkgModelPass, 0, sigTemplates, 0)

                    # ROOT.set_massRange(MassMinSta_, MassMaxSta_, MassBinSta_)
                    # ROOT.getZyield(hGlopass, m, "Glo", etaRegion, sigModel, bkgModelPass, 1, sigTemplates, 0)
                    # ROOT.getZyield(hGlofail, m, "Glo", etaRegion, sigModel, bkgModelFail, 0, sigTemplates, 0)
                    # ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

                    # ROOT.getZyield(hStapass, m, "Sta", etaRegion, sigModel, bkgModelPass, 1, sigTemplates, 0)
                    # ROOT.getZyield(hStafail, m, "Sta", etaRegion, sigModel, bkgModelFail, 0, sigTemplates, 0)

                    # remove the histogram templates, not needed anymore
                    os.system("rm {0}/histTemplates_*".format(outSubDir))

            # cHLT = ROOT.extractCorrelation_HLT(sigTemplates, hPV, etaRegionZ)
            # result = extract_results(outSubDir, m, 1, hPV, mcCorrelations)

            loaded_histograms = {"HLT2": h2HLT, "HLT1": h1HLT, "IDfail": hIDfail, "Glopass": hGlopass, "Glofail": hGlofail, "Stapass": hStapass, "Stafail": hStafail}
            # massbins = {"Other_Categories": (MassBin_, MassMin_, MassMax_), "Glo": (MassBinSta_, MassMinSta_, MassMaxSta_)}

            
            # load templates, correlations, acceptance
            # load reco-cut histograms
            templates_reco_hdf5 = h5py.File('templates/templates_reco.hdf5', 'r')

            templates_reco = {}
            for key in templates_reco_hdf5.keys():
                templates_reco[key] = np.array(templates_reco_hdf5[key])

            # reduce to 1d template in mass
            nPV = hPV.view()
            templates_1d = {}

            for key in templates_reco.keys():
                templates_1d[key] = np.matmul(templates_reco[key], nPV)

            # for key in templates_reco.keys():
            #     templates_1d[key] = np.sum(templates_reco[key], axis = 1)
            
  
              # load correlation templates
            correlation_hdf5 = h5py.File("templates/correlations.hdf5", 'r')

            correlation_templates = {}
            for key in correlation_hdf5.keys():
                correlation_templates[key] = np.array(correlation_hdf5[key])
                correlation_templates[key] = np.nan_to_num(correlation_templates[key], nan=1)


            # load acceptance templates
            acceptance_hdf5 = h5py.File("templates/acceptance.hdf5", 'r')

            acceptance_templates = {}
            for key in acceptance_hdf5.keys():
                acceptance_templates[key] = np.array(acceptance_hdf5[key])

            pdb.set_trace()

            result = extract_results(loaded_histograms, templates_1d, correlation_templates, acceptance_templates, hPV)
            

            pdb.set_trace()

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

                # delivered, efficiency corrected Z boson count extrapolated using reference luminosity
                delZCount = result["recZCount"] / recLumi * delLumi

                # instantaneous delivered Z rate (corrected for deadtime)
                ZRate = result["recZCount"] / (timewindow * deadtime)

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
                    "instDelLumi": delLumi/timewindow,
                    "ZRate": ZRate
                })
            
                results.append(result)
            else:
                log.info("No result - continue")
        
        if measurement is None or measurement == m:
            utils.writePerRunCSV(results, outCSVDir, run)
            

        results = []


    #if args.writeSummaryCSV:
    #    utils.writeSummaryCSV(outCSVDir)
    #    utils.writeSummaryCSV(outCSVDir, outName="cms_2.4", keys=["fill","run"])
    
    # To build the merge csv file locally 
    utils.writeSummaryCSV(outCSVDir)
    utils.writeSummaryCSV(outCSVDir, outName="cms_2.4", keys=['fill','beginTime','endTime','ZRate','instDelLumi','delLumi','delZCount'])
    print("Done!") # needed for submit script to check if job is completed
