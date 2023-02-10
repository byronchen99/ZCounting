import ROOT
import pandas as pd
import glob
import os
import pdb
import uncertainties as unc
import datetime
import numpy as np

from python import utils, logging, common

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None

# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

# -------------------------------------------------------------------------------------
# helper function to read the workspace for a specific fit
def open_workspace(directory, filename, m):
    file_result = f"{directory}/{filename}_{m}.root"
    
    if not os.path.isfile(file_result):
        log.warning("No result for `{0}`".format(file_result))
        return None, None
    
    f = ROOT.TFile(file_result,"READ")
    w = f.Get("workspace")
    return f, w

def open_workspace_yield(directory, filename, m):
    f, w = open_workspace(directory, "workspace_yield_"+filename, m)
    
    Nsig = w.var("Nsig").getVal()
    chi2 = w.arg("chi2").getVal()
    
    f.Close()
    f.Delete()
    w.Delete()
    return Nsig, chi2

def unorm(value):
    if value > 0:
        return unc.ufloat(value, np.sqrt(value))
    else: 
        return unc.ufloat(value, value)

# -------------------------------------------------------------------------------------
def extract_results(directory, m, cIO, cID, cHLT, cKinematicSelection):
    log.info("Extracting fit results in {0} for {1}".format(directory,m))  

    # --- For identification (ID) efficiency        
    NsigHLT2, chi2HLT2 = open_workspace_yield(directory, "HLT_{0}_2".format(etaRegion), m)
    NsigHLT1, chi2HLT1 = open_workspace_yield(directory, "HLT_{0}_1".format(etaRegion), m)
    NsigIDFail, chi2ID = open_workspace_yield(directory, "ID_{0}_0".format(etaRegion), m)

    NsigHLT2 = unorm(NsigHLT2)
    NsigHLT1 = unorm(NsigHLT1)
    NsigIDFail = unorm(NsigIDFail)

    effHLT = (2 * NsigHLT2) / (2 * NsigHLT2 + NsigHLT1)
    effID = (2 * NsigHLT2 + NsigHLT1) / (2 * NsigHLT2 + NsigHLT1 + NsigIDFail)

    NsigGloPass, chi2GloPass = open_workspace_yield(directory, "Glo_{0}_1".format(etaRegion), m)
    NsigGloFail, chi2GloFail = open_workspace_yield(directory, "Glo_{0}_0".format(etaRegion), m)

    NsigStaPass, chi2StaPass = open_workspace_yield(directory, "Sta_{0}_1".format(etaRegion), m)
    NsigStaFail, chi2StaFail = open_workspace_yield(directory, "Sta_{0}_0".format(etaRegion), m)

    NsigGloPass = unorm(NsigGloPass)
    NsigGloFail = unorm(NsigGloFail)

    NsigStaPass = unorm(NsigStaPass)
    NsigStaFail = unorm(NsigStaFail)

    effGlo = NsigGloPass / (NsigGloPass + NsigGloFail)
    effSta = NsigStaPass / (NsigStaPass + NsigStaFail)

    effMu = effID*effGlo*effSta

    res = {
        "NsigHLT2": NsigHLT2,
        "NsigHLT1": NsigHLT1,
        "NsigIDFail": NsigIDFail,
        "NsigGloPass": NsigGloPass,
        "NsigGloFail": NsigGloFail,
        "NsigStaPass": NsigStaPass,
        "NsigStaFail": NsigStaFail,
        "chi2HLT2": chi2HLT2,
        "chi2HLT1": chi2HLT1,
        "chi2GloPass": chi2GloPass,
        "chi2GloFail": chi2GloFail,
        "chi2StaPass": chi2StaPass,
        "chi2StaFail": chi2StaFail,
        "effHLT": effHLT,
        "effID": effID,
        "effGlo": effGlo,
        "effSta": effSta,
        "effMu": effMu,
        "cHLT": cHLT,
        "cIO": cIO,
        "cID": cID,
        "cKinematicSelection": cKinematicSelection,
        "recZCount": (NsigHLT2 + 0.5*NsigHLT1)**2/NsigHLT2 / effMu**2 * cHLT * cID * cIO**2 * cKinematicSelection
    }

    return res

# -------------------------------------------------------------------------------------
def extract_results_dqm(directory, m, cHLT):
    log.info("Extracting fit results in {0} for {1}".format(directory,m))  

    # --- For identification (ID) efficiency        
    NsigHLT2, chi2HLT2 = open_workspace_yield(directory, "HLT_{0}_2".format(etaRegion), m)
    NsigHLT1, chi2HLT1 = open_workspace_yield(directory, "HLT_{0}_1".format(etaRegion), m)
    NsigIDFail, chi2ID = open_workspace_yield(directory, "Sel_{0}_0".format(etaRegion), m)

    NsigHLT2 = unorm(NsigHLT2)
    NsigHLT1 = unorm(NsigHLT1)
    NsigIDFail = unorm(NsigIDFail)

    effHLT = (2 * NsigHLT2) / (2 * NsigHLT2 + NsigHLT1)
    effID = (2 * NsigHLT2 + NsigHLT1) / (2 * NsigHLT2 + NsigHLT1 + NsigIDFail)

    NsigGloFail, chi2GloFail = open_workspace_yield(directory, "Glo_{0}_0".format(etaRegion), m)

    NsigGloFail = unorm(NsigGloFail)

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
    import os

    cmsswbase = os.environ['CMSSW_BASE']
    # cmsswbase = "/eos/home-d/dwalter/Lumi/ZCounting_Standalone/CMSSW_12_4_10/src/ZCounting/ZHarvester"

    parser = common.parser()
    parser = common.parser_zharvest(parser)

    parser.add_argument("-m", "--measurement", type=int, default=None, 
        help="Only fit a specific measurement of the run")
    args = parser.parse_args()

    log = logging.setup_logger(__file__, args.verbose)

    prefix_dqm="ZCountingInOut-V17_38-"

    ########################################
    # link to resouces
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
        byLsCSV = f"{args.input}/2022/brilcalcByLS/byLS_Collisions22_355100_362760_Muon_20230210.csv"
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

    if args.mode == "DQM":
        npvMin_ = 0.5     
        npvMax_ = 100.5   
    else:
        npvMin_ = -0.5
        npvMax_ = 74.5

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

    if args.mode == "TTrees":
        # Define histograms for fitting
        hPV = ROOT.TH1D("h_PV","", npvBin_, npvMin_, npvMax_)

        h2HLT = ROOT.TH1D("h_mass_2HLT_Z","",MassBin_, MassMin_, MassMax_)
        h1HLT = ROOT.TH1D("h_mass_1HLT_Z","",MassBin_, MassMin_, MassMax_)
        hIDfail = ROOT.TH1D("h_mass_ID_fail","",MassBin_, MassMin_, MassMax_)

        hGlopass = ROOT.TH1D("h_mass_Glo_pass","",MassBinSta_, MassMinSta_, MassMaxSta_)
        hGlofail = ROOT.TH1D("h_mass_Glo_fail","",MassBinSta_, MassMinSta_, MassMaxSta_)

        hStapass = ROOT.TH1D("h_mass_Sta_pass","",MassBin_, MassMin_, MassMax_)
        hStafail = ROOT.TH1D("h_mass_Sta_fail","",MassBin_, MassMin_, MassMax_)
    
    #####################################   
    recLumi = 0
    firstRun = 0
    lastRun = 0
    df=None
    results = []
    mergeNextRun=False
    log.info(f"Looping over runs... {args.beginRun} to {args.endRun}")
    for run, byLS_run in byLS_data.groupby('run', sort=True):
        
        # first and last run of the measurement
        if firstRun == 0:
            firstRun = run
        lastRun = run

        fill = byLS_run.drop_duplicates('fill')['fill'].values[0]
        LSlist = byLS_run['ls'].values.tolist()
        Lumilist = byLS_run.loc[byLS_run['ls'].isin(LSlist)]['recorded(/pb)'].values.tolist()

        log.info(f"Running Fill {fill}")
        log.info(f"Running Run {run}")
        
        if args.mode == "TTrees":
            eosFile = f"{eosDir}/{prefix_dqm}*Muon_{run}*.root"
            eosFiles = glob.glob(eosFile)
            if len(eosFiles) == 1:
                eosFile = eosFiles[0]
            else:
                log.warning("No file or more than one was found! - continue")
                log.warning(f"Was looking for: {eosFile}")            
                continue
            file_ = ROOT.TFile(eosFile,"READ")

            # histograms need to be in same directory so that they can get filled
            hPV.SetDirectory(file_)
            h2HLT.SetDirectory(file_)
            h1HLT.SetDirectory(file_)
            hIDfail.SetDirectory(file_)
            hGlopass.SetDirectory(file_)
            hGlofail.SetDirectory(file_)
            hStapass.SetDirectory(file_)
            hStafail.SetDirectory(file_)

            # trees with muon pairs
            tHLT = file_.Get("HLT")
            tID  = file_.Get("ID")
            tGlo = file_.Get("Glo")
            tSta = file_.Get("Sta")

            # ZCountlist = [tHLT.GetEntries("lumiBlock=={0}".format(l)) for l in LSlist]
        else:
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
            log.debug("Selected lumi section list {0}".format(goodLSlist))

            if measurement is not None and measurement < m:
                break
                    
            # create datafram byLS for measurement
            byLS_m = byLS_run.loc[byLS_run['ls'].isin(goodLSlist)]
            
            if args.mode == "TTrees":
                ### fill histograms
                file_.cd() # switch to directory where ttrees and histograms are placed

                # define acceptance cuts
                acceptance = " && mass>={0} && mass<{1} && ptTag > {2} && ptProbe > {2} && abs(etaTag) < {3} && abs(etaProbe) < {3}".format(MassMin_, MassMax_, args.ptCut, args.etaCut)
                acceptanceSta = " && mass>={0} && mass<{1} && ptTag > {2} && ptProbe > {2} && abs(etaTag) < {3} && abs(etaProbe) < {3}".format(MassMinSta_, MassMaxSta_, args.ptCut, args.etaCut)

                log.info("Fill histograms for measurement {0} ...".format(m))                        
                for iLS in goodLSlist:
                        
                    tHLT.Draw("nPV>>+h_PV","lumiBlock=={0}".format(iLS))
                    
                    n2Before = h2HLT.Integral()
                    n1Before = h1HLT.Integral()

                    tHLT.Draw("mass>>+h_mass_2HLT_Z", "pass==2 && lumiBlock=={0} {1}".format(iLS, acceptance))
                    tHLT.Draw("mass>>+h_mass_1HLT_Z", "pass==1 && lumiBlock=={0} {1}".format(iLS, acceptance))
                    tID.Draw("mass>>+h_mass_ID_fail", "pass==0 && lumiBlock=={0} {1}".format(iLS, acceptance))
                    
                    tGlo.Draw("mass>>+h_mass_Glo_pass","pass==1 && lumiBlock=={0} {1}".format(iLS, acceptanceSta))                
                    tGlo.Draw("mass>>+h_mass_Glo_fail","pass==0 && lumiBlock=={0} {1}".format(iLS, acceptanceSta))                

                    tSta.Draw("mass>>+h_mass_Sta_pass","pass==1 && lumiBlock=={0} {1}".format(iLS, acceptance))   
                    tSta.Draw("mass>>+h_mass_Sta_fail","pass==0 && lumiBlock=={0} {1}".format(iLS, acceptance))  

                    n2After = h2HLT.Integral()
                    n1After = h1HLT.Integral()
                    
                    # store the number of 1hlt and 2hlt events in each lumisection
                    n2 = n2After - n2Before
                    n1 = n1After - n1Before

                    byLS_m.loc[byLS_m['ls'] == iLS, 'N2HLT'] = n2
                    byLS_m.loc[byLS_m['ls'] == iLS, 'N1HLT'] = n1                    
            else:
                # get histograms from good lumisections
                log.info("Load histograms ...")
                if args.mode == "DQM":
                    prefix=f"DQMData/Run {run}/ZCounting/Run summary/Histograms/"
                else:
                    prefix=""
                
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

            # mask_HLT2 = "pass==2 && "+mask_lumi+acceptance
            # mask_HLT1 = "pass==1 && "+mask_lumi+acceptance

            # dHLT2 = dHLT.Filter("pass==2 && "+mask_lumi+acceptance).AsNumpy(["mass"])
            # dHLT1 = dHLT.Filter("pass==1 && "+mask_lumi+acceptance).AsNumpy(["mass"])

            # mass = ROOT.RooRealVar("mass", "mass", 60, 120)

            # data = ROOT.RooDataSet.from_numpy({"x": dHLT2}, [mass])

            # read ttree into RooDataSet
            # dHLT2 = dHLT.Filter("pass==2 && "+mask_lumi+acceptance).Take['float']("mass")
            # dHLT1 = dHLT.Filter("pass==2 && "+mask_lumi+acceptance).Take['float']("mass")


            if df is None:
                df = byLS_m
            else:
                df = df.append(byLS_m, sort=False)
            
            recLumi = df['recorded(/pb)'].sum()

            log.info("Have now recorded lumi = {0}".format(recLumi))            
            # log.info("Have now {0} | {1} events".format(df['N2HLT'].sum(), h2HLT.Integral()))
            
            if args.mode == "TTrees":
                # check if upcoming runs make enough data for a measurement
                lumi=0
                mergeNextRun=True
                nextRun = run+1
                for r, dr in byLS_data.groupby('run'):
                    if r <= run or r >= int(args.endRun):
                        continue
                    # if utils.getFileName(eosDir, r) is None: # check if file of next run exists
                    #     continue
                    nextRun = r
                    LS = dr['ls'].values.tolist()
                    lumi += sum(dr.loc[dr['ls'].isin(LS)]['recorded(/pb)'].values)
                    if lumi > 0.5 * LumiPerMeasurement:
                        mergeNextRun = False
                        break
                
                mergeNextRun = nextRun < int(args.endRun) and (mergeNextRun or recLumi < 0.5 * LumiPerMeasurement)            

                if mergeNextRun:
                    log.info("Merge with next run ... ")
                    continue

            log.info("Histograms filled ...")  

            if firstRun != lastRun:
                outSubDir = outDir + "Run{0}to{1}".format(firstRun,lastRun)
            else:
                outSubDir = outDir + "Run{0}/".format(run)

            log.debug("Running measurement {0}".format(m))

            ### --- perform fit with zfit
            # from python.fit import zfit

            # pdb.set_trace()

            # h2HLT = utils.np_to_hist(h2HLT, MassBin_, MassMin_, MassMax_, "hist")

            # fit(h2HLT, MassBin_, MassMin_, MassMax_, category="HLT 2")
            # fit(h1HLT, MassBin_, MassMin_, MassMax_, category="HLT 1")
            # fit(hIDfail, MassBin_, MassMin_, MassMax_, category="ID fail")
            # fit(hGlofail, MassBin_, MassMinSta_, MassMaxSta_, category="Glo fail")

            ### --- 

            if not args.collect:
                
                if measurement is None or measurement == m:
                    # skip the fit if we look for another measurement
                
                    if not os.path.isdir(outSubDir):
                        os.mkdir(outSubDir)
                    
                    ROOT.set_output(outSubDir)
                    ROOT.set_luminosity(recLumi)

                    if args.mode == "TTrees":
                        ROOT.getZyield(h2HLT, m, "HLT", etaRegion, sigModel, bkgModelPass, 2, sigTemplates, 0)
                        ROOT.getZyield(h1HLT, m, "HLT", etaRegion, sigModel, bkgModelPass, 1, sigTemplates, 0)
                        ROOT.getZyield(hIDfail, m, "ID", etaRegion, sigModel, bkgModelPass, 0, sigTemplates, 0)

                        ROOT.set_massRange(MassMinSta_, MassMaxSta_, MassBinSta_)
                        ROOT.getZyield(hGlopass, m, "Glo", etaRegion, sigModel, bkgModelPass, 1, sigTemplates, 0)
                        ROOT.getZyield(hGlofail, m, "Glo", etaRegion, sigModel, bkgModelFail, 0, sigTemplates, 0)
                        ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

                        ROOT.getZyield(hStapass, m, "Sta", etaRegion, sigModel, bkgModelPass, 1, sigTemplates, 0)
                        ROOT.getZyield(hStafail, m, "Sta", etaRegion, sigModel, bkgModelFail, 0, sigTemplates, 0)                        

                    else:
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

            if args.mode == "TTrees":
                cID = utils.getCorrelation(hPV, mcCorrelations, which="ID")
                cIO = utils.getCorrelation(hPV, mcCorrelations, which="IO")
                cKinematicSelection = utils.getCorrelation(hPV, mcCorrelations, which="KinematicSelection")
                result = extract_results(outSubDir, m, cIO, cID, cHLT, cKinematicSelection)
            else:
                result = extract_results_dqm(outSubDir, m, cHLT)

            if result:              
                delLumi = df['delivered(/pb)'].sum()
                recLumi = df['recorded(/pb)'].sum()
                # deadtime - fraction during the good lumisections
                deadtime = recLumi / delLumi
                timewindow = len(df) * secPerLS

                # convert time string to datetime format
                beginTime = utils.to_DateTime(df['time'][0], string_format = "mm/dd/yy")
                endTime = utils.to_DateTime(df['time'][-1], string_format = "mm/dd/yy") + datetime.timedelta(seconds=secPerLS)
                
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
                    "pileUp": df['avgpu'].mean(),
                    "delZCount": delZCount,
                })
            
                results.append(result)
            else:
                log.info("No result - continue")
            
            # prepare for next measurement
            df=None

            if args.mode == "TTrees":
                # clean the histograms for the next measurement
                hPV.Reset() 
                h2HLT.Reset()
                h1HLT.Reset()
                hIDfail.Reset()

                hGlopass.Reset()
                hGlofail.Reset()
                hStapass.Reset()
                hStafail.Reset()


        ### prepare for next run
        if args.mode == "TTrees":
            # keep histograms
            hPV.SetDirectory(0)
            h2HLT.SetDirectory(0)
            h1HLT.SetDirectory(0)
            hIDfail.SetDirectory(0)

            hGlopass.SetDirectory(0)
            hGlofail.SetDirectory(0)
            hStapass.SetDirectory(0)
            hStafail.SetDirectory(0)

            file_.Close()

        if mergeNextRun:
            continue
        
        if measurement is None or measurement == m:
            ## Write per measurement csv file - one per run
            log.info("Writing per Run CSV file")
            results = pd.concat([pd.DataFrame([result]) for result in results], ignore_index=True, sort=False)

            with open(outCSVDir + '/csvfile{0}.csv'.format(run), 'w') as file:
                results.to_csv(file, index=False)

        firstRun = 0
        results = []

    if args.writeSummaryCSV:
        utils.writeSummaryCSV(outCSVDir, writeByLS=False)

    log.info("Done")
