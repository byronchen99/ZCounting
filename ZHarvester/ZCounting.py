import logging as log
import ROOT
import pandas as pd
import glob
import os
import pdb
import uncertainties as unc
import datetime
import numpy as np

from python.utils import writeSummaryCSV, getEra, getFileName, load_input_csv, get_ls_for_next_measurement, getCorrelation, to_DateTime


# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None

# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

def extract_results(directory, m, cIO, cID):
    log.info(" === Extracting fit results in {0} for {1}".format(directory,m))
    
    # helper function to read the workspace for a specific fit
    def open_workspace(filename):
        file_result = directory+"/{0}_{1}.root".format(filename, m)
        
        if not os.path.isfile(file_result):
            log.warning(" === No result for `{0}`".format(file_result))
            return None, None
        
        f = ROOT.TFile(file_result,"READ")
        w = f.Get("workspace")
        return f, w

    def unorm(value):
        if value > 0:
            return unc.ufloat(value, np.sqrt(value))
        else: 
            return unc.ufloat(value, value)

    # --- HLT efficiency and yield
    f, w = open_workspace("workspace_{0}".format(etaRegionZ))
    
    effHLT = w.var("eff").getVal()
    cHLT = w.arg("c").getVal()
    NsigHLT = w.var("Nsig").getVal()
    chi2HLT = w.arg("chi2").getVal()
    
    f.Close()
    f.Delete()
    w.Delete()

    NsigHLT1 = unorm(2*effHLT*(1-cHLT*effHLT) * NsigHLT)
    NsigHLT2 = unorm(effHLT**2 * cHLT * NsigHLT)
    effHLT = (2 * NsigHLT2) / (2 * NsigHLT2 + NsigHLT1)

    res = {
            "effHLT": effHLT,
            "NsigHLT1": NsigHLT1,
            "NsigHLT2": NsigHLT2,
            "cHLT": cHLT,
            "cIO": cIO,
            "chi2HLT": chi2HLT,
            "cID": cID,
    }

    if factorization_scheme == 1:
        # --- For selection (Sel) efficiency        
        f, w = open_workspace("workspace_Sel_{0}".format(etaRegion))
        
        effSel = w.var("eff").getVal()
        NsigSel = w.var("Nsig").getVal()
        chi2Sel = w.arg("chi2").getVal() 

        f.Close()
        f.Delete()
        w.Delete()
    
        # --- For global (Glo) efficiency
        f, w = open_workspace("workspace_Glo_{0}".format(etaRegion))
     
        effGlo = w.var("eff").getVal()
        NsigGlo = w.var("Nsig").getVal()
        chi2Glo = w.arg("chi2").getVal()

        f.Close()
        f.Delete()
        w.Delete()

        # --- For standalone (Sta) efficiency
        f, w = open_workspace("workspace_Sta_{0}".format(etaRegion))

        effSta = w.var("eff").getVal()
        NsigSta = w.var("Nsig").getVal()
        chi2Sta = w.arg("chi2").getVal()

        f.Close()
        f.Delete()
        w.Delete()

        # --- Final calculation, poisson uncertainty from counts 
        NsigSelPass = unorm(NsigSel*effSel)
        NsigSelFail = unorm(NsigSel*(1-effSel))

        NsigGloPass = unorm(NsigGlo*effGlo)
        NsigGloFail = unorm(NsigGlo*(1-effGlo))

        NsigStaPass = unorm(NsigSta*effSta)
        NsigStaFail = unorm(NsigSta*(1-effSta))

        effSel = NsigSelPass / (NsigSelPass + NsigSelFail)
        effGlo = NsigGloPass / (NsigGloPass + NsigGloFail)
        effSta = NsigStaPass / (NsigStaPass + NsigStaFail)

        effID = (effSta*effGlo*effSel)**2

        res.update({
            "effSel": effSel,
            "effGlo": effGlo,
            "effSta": effSta,
            "NsigSelPass": NsigSelPass,
            "NsigSelFail": NsigSelFail,
            "NsigGloPass": NsigGloPass,
            "NsigGloFail": NsigGloFail,
            "NsigStaPass": NsigStaPass,
            "NsigStaFail": NsigStaFail,
            "chi2Sel": chi2Sel,
            "chi2Glo": chi2Glo,
            "chi2Sta": chi2Sta,
        })

    else:
        # --- For selection (Sel) efficiency        
        f, w = open_workspace("workspace_yield_Sel_{0}_0".format(etaRegion))
        
        NsigSelFail = w.var("Nsig").getVal()
        chi2Sel = w.arg("chi2").getVal()
        
        f.Close()
        f.Delete()
        w.Delete()

        # --- For tracking (Trk) efficiency        
        f, w = open_workspace("workspace_Trk_{0}".format(etaRegion))
        
        effTrk = w.var("eff").getVal()
        NsigTrk = w.var("Nsig").getVal()
        chi2Trk = w.arg("chi2").getVal()

        f.Close()
        f.Delete()

        # # --- For corrections of (Trk) efficiency
        # f, w = open_workspace("workspace_yield_Trk_{0}_2".format(etaRegion))

        # # in this case it is possible that there was no event and the fit was not performed    
        # if w is None:
        #     NsigTrkPass2 = 0
        #     chi2TrkPass2 = -1 
        # else:
        #     NsigTrkPass2 = w.var("Nsig").getVal()
        #     chi2TrkPass2 = w.arg("chi2").getVal() 
            
        #     f.Close()
        #     f.Delete()
        #     w.Delete()

        # --- Final calculation, poisson uncertainty from counts 
        NsigSelFail = unorm(NsigSelFail)

        NsigTrkPass = unorm(NsigTrk*effTrk)
        NsigTrkFail = unorm(NsigTrk*(1-effTrk))

        # NsigTrkPass2 = unorm(NsigTrkPass2)

        effSel = (2 * NsigHLT2 + NsigHLT1) / (2 * NsigHLT2 + NsigHLT1 + NsigSelFail)

        # # tracking efficiency: corrections from fakes
        # # first order correction
        # corr1 = NsigTrkPass2 * NsigTrkFail / NsigTrkPass
        # # second order correction
        # corr2 = NsigTrkPass2**2 * NsigTrkFail**2 / NsigTrkPass**3
        # effTrk = (NsigTrkPass + NsigTrkPass2 - corr1 + corr2) / (NsigTrkPass + NsigTrkPass2 + NsigTrkFail)

        effTrk = NsigTrkPass / (NsigTrkPass + NsigTrkFail)


        res.update({
            "effSel": effSel,
            "effTrk": effTrk,
            "NsigSelFail": NsigSelFail,
            # "NsigTrkPass2": NsigTrkPass2,
            "NsigTrkPass": NsigTrkPass,
            "NsigTrkFail": NsigTrkFail,
            "chi2Sel": chi2Sel,
            "chi2Trk": chi2Trk,
            # "chi2TrkPass2": chi2TrkPass2,
        })

        effID = (effSel*effTrk)**2

    res.update({
        "recZCount": (NsigHLT2 + 0.5*NsigHLT1)**2/NsigHLT2 * cHLT * cID * cIO**2 / effID
    })

    return res


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
    parser.add_argument("-v", "--verbose", help="increase logging level from INFO to DEBUG", default=False,
                        action="store_true")
    parser.add_argument("-c", "--writeSummaryCSV", default=False, action="store_true",
                        help="produce merged CSV with all runs")
    parser.add_argument("-i", "--dirDQM", help="Directory to the input root files from the DQM Offline module",
                        default="default")
    parser.add_argument("--byLsCSV", help="ByLs csv input generated by testBril.sh",
                        default="default")
    parser.add_argument("--sigTemplates", default="default", type=str,
        help="Choose one of the options for signal model (MC, MCxGaus, MCxCB, BW, BWxCB, BWxGaus). Default is MCxGaus and exp")
    parser.add_argument("--bkgTemplates", default="default", type=str,
        help="Choose one of the options for background model (Exp, Quad, QuadPlusExp, CMSShape, Das). Default is CMSShape and linear")
    parser.add_argument('--ptCut', type=float, default=25.,
                        help='specify lower pt cut on tag and probe muons')
    parser.add_argument('--etaCut', type=float, default=2.4,
                        help='specify upper |eta| cut on tag and probe muons')
    parser.add_argument('--mass', nargs=3, metavar=('LOW', 'HIGH', 'NUMBER'), default=(60,120,120), type=int,
                        help='specify mass range for tag and probe muon pairs')
    parser.add_argument('--LumiPerMeasurement', default=20, type=float,
                        help='specify amount of luminosity per measurement in pb-1')
    parser.add_argument('--inclusive', default=False, action="store_true",
                        help='specify whether or not to do an inclusive fit of the specified runs')
    parser.add_argument('--collect', default=False, action="store_true",
                        help='specify whether or not to run the fits or just collect the results')
    parser.add_argument('-f','--factorization', default=2, type=int,
                        help='specify factorization schemes for muon efficiency')
    parser.add_argument("-o", "--dirOut", help="where to store the output files", default="./")

    args = parser.parse_args()
    if args.verbose:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
    else:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.INFO)

    # factorization schemes for efficiency of muon
    # implemented so far:
    #    1: eff(Sel|Glo) * eff(Glo|Sta) * eff(Sta|Trk)
    #    2: eff(Sel|Trk) * eff(Trk|Sta)
    factorization_scheme = args.factorization

    if factorization_scheme == 1:
        mcCorrelations = "/mcCorrections/InnerOuter_V17_01/"
        prefix_dqm="ZCountingAll-V17_02-"
    elif factorization_scheme == 2:
        mcCorrelations = "/mcCorrections/correlations_V17_31/"
        prefix_dqm="ZCountingAll-V17_31-"

    ########################################
    # link to resouces
    eosDir           = args.dirDQM
     #"DQMData/Run {0}/ZCounting/Run summary/Histograms/".format(run)
    resPath = cmsswbase + "/src/ZCounting/ZHarvester/res/"
    if( args.beginRun >= 272007 and args.beginRun < 278808
        # there is an overlap for 2016 F in runs with pre and post VFP settings
        and args.beginRun not in [278769, 278801, 278802, 278803, 278804, 278805, 278808]
    ):                                                          # 2016 pre VFP
        currentYear = 2016
        byLsCSV          = resPath+"/FillByLs_2016.csv"
        mcCorrelations   = resPath+mcCorrelations+"c_nPV_2016preVFP.root"
        sigTemplates     = eosDir+"/"+prefix_dqm+"Summer16preVFP-DYJetsToLL_M_50_LO.root"
    elif args.beginRun < 294645:                                # 2016 post VFP
        currentYear = 2016
        byLsCSV          = resPath+"/FillByLs_2016.csv"
        mcCorrelations   = resPath+mcCorrelations+"c_nPV_2016postVFP.root"
        sigTemplates     = eosDir+"/"+prefix_dqm+"Summer16postVFP-DYJetsToLL_M_50_LO.root"
    elif args.beginRun > 297020 and args.beginRun < 306828:     # 2017
        currentYear = 2017
        byLsCSV          = resPath+"/FillByLs_2017_IsoMu24.csv"
        mcCorrelations   = resPath+mcCorrelations+"c_nPV_2017.root"
        sigTemplates     = eosDir+"/"+prefix_dqm+"Fall17-DYJetsToLL_M_50_LO.root"
    elif args.beginRun >= 306926 and args.beginRun < 307083:    # 2017 H
        currentYear = 2017
        byLsCSV          = resPath+"/FillByLs_2017_lowPU.csv"
        mcCorrelations   = resPath+mcCorrelations+"c_nPV_2017.root"
        sigTemplates     = eosDir+"/"+prefix_dqm+"Fall17-DYJetsToLL_M_50_LO.root"
    elif args.beginRun >= 315252 and args.beginRun < 325273:    # 2018
        currentYear = 2018
        byLsCSV          = resPath+"/FillByLs_2018.csv"
        mcCorrelations   = resPath+mcCorrelations+"c_nPV_2018.root"
        sigTemplates     = eosDir+"/"+prefix_dqm+"Autumn18-DYJetsToLL_M_50_LO.root"
    elif args.beginRun >= 355100:                               # 2022
        currentYear = 2022
        byLsCSV = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/brilcalcByLS/byLS_Collisions22_355100_356615_Golden.csv"
        mcCorrelations   = "/eos/cms/store//group/comm_luminosity/ZCounting/2022/CorrelationFactors/c_nPV_2022.root"
        prefix_dqm =  "ZCountingAll-V01-"
        sigTemplates = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/SignalTemplates/ZCountingAll-V01-Winter22-DYJetsToLL_M_50_LO.root"
    else:
        currentYear = 2017
        mcCorrelations      = None
        byLsCSV             = None
        sigTemplates        = None

    byLsCSV          = byLsCSV          if args.byLsCSV       == "default"   else args.byLsCSV
    measurement      = args.measurement

    log.info("----------------------------------")
    log.info("Use eosDir:              {0}".format(eosDir))
    log.info("Use byLsCSV:             {0}".format(byLsCSV))
    log.info("Use sigTemplates:        {0}".format(sigTemplates))
    log.info("Mass range from:         {0} to {1}".format(*args.mass))
    log.info("Lumi per Measurement:    {0}".format(args.LumiPerMeasurement))
    log.info("----------------------------------")
    
    # signal model
    if args.sigTemplates == "MCxGauss" or args.sigTemplates == "default":
        sigModel = 2 # MC, folding with gauss
    elif args.sigTemplates == "MC":
        sigModel = 4 # MC, no folding
    elif args.sigTemplates == "BW":
        sigModel = 3 # BW, no folding
    elif args.sigTemplates == "BWxCB":
        sigModel = 1 # BW, folding with crystal ball
    elif args.sigTemplates == "BWxGaus":
        sigModel = 5 # BW, folding with gauss
    elif args.sigTemplates == "MCxCB":
        sigModel = 6 # MC, folding with crystal ball
    else:
        log.warning("signal model {0} unknown! exit()".format(args.sigTemplates))
        exit()

    # background model 
    if args.bkgTemplates == "default" :
        bkgModelPass = 1    # exp
        bkgModelFail = 6    # RooCMSShape
    elif args.bkgTemplates == "alt" :
        bkgModelPass = 8    # linear
        bkgModelFail = 4    # Das
    elif args.bkgTemplates == "Exp":
        bkgModelFail = 1
        bkgModelPass = 1
    elif args.bkgTemplates == "Quad":
        bkgModelFail = 2
        bkgModelPass = 2
    elif args.bkgTemplates == "QuadPlusExp":
        bkgModelFail = 3
        bkgModelPass = 3
    elif args.bkgTemplates == "Das":
        bkgModelFail = 4
        bkgModelPass = 4
    elif args.bkgTemplates == "CMSShape":
        bkgModelFail = 6
        bkgModelPass = 6
    else:
        log.warning("background model {0} unknown! exit()".format(args.bkgTemplates))
        exit()
        
    ########## Input configuration ##########
    # ByLS csv inputs generated by testBRIL.sh
    byLS_filelist = glob.glob(byLsCSV)
    byLS_filelist.sort(key=os.path.getmtime)
    byLS_filename = byLS_filelist[-1]
    log.info(" The brilcalc csv file: " + str(byLS_filename))


    outDir = args.dirOut if args.dirOut.endswith("/") else args.dirOut+"/"
    if not os.path.isdir(outDir):
        os.mkdir(outDir)

    outCSVDir = outDir+"csvFiles/"
    if not os.path.isdir(outCSVDir):
        try:
            os.mkdir(outCSVDir)
        except OSError:
            log.warning(": directory already exists ...")

    ########### Constant settings ##########
    secPerLS = float(23.3)
    
    LumiPerMeasurement = args.LumiPerMeasurement  # minimum recorded lumi for one measurement in pb-1

    #configuration of fit models
    MassMin_ = int(args.mass[0])
    MassMax_ = int(args.mass[1])
    MassBin_ = int(args.mass[2])

    MassMinSta_ = int(50)
    MassMaxSta_ = int(130)
    MassBinSta_ = int(MassMaxSta_ - MassMinSta_)*2

    npvBin_ = 75
    npvMin_ = -0.5
    npvMax_ = 74.5

    if args.etaCut == 0.9:
        etaRegion = "B"
        etaRegionZ = "BB"
    else:
        etaRegion = "I"
        etaRegionZ = "I"

    if not args.collect:
        log.info(" Loading C marco...")
        # load functions for fitting
        ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(
            __file__)) + "/calculateDataEfficiency.C")

        ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

        ROOT.set_npvRange(npvMin_, npvMax_)
        if currentYear >= 2022:
            ROOT.set_energy(13.6)

        ROOT.set_ptCut(args.ptCut)
        ROOT.set_etaCut(args.etaCut)

        # ROOT.set_IntegrateBinsPrecision(0.1)
    
    byLS_data = load_input_csv(byLS_filename)

    #####################################   

    # For fitting
    hPV = ROOT.TH1D("h_PV","", npvBin_, npvMin_, npvMax_)

    h2HLT = ROOT.TH1D("h_mass_2HLT_Z","",MassBin_, MassMin_, MassMax_)
    h1HLT = ROOT.TH1D("h_mass_1HLT_Z","",MassBin_, MassMin_, MassMax_)
    hSelfail = ROOT.TH1D("h_mass_Sel_fail","",MassBin_, MassMin_, MassMax_)

    if factorization_scheme == 1:
        hSelpass = ROOT.TH1D("h_mass_Sel_pass","",MassBin_, MassMin_, MassMax_)

        hGlopass = ROOT.TH1D("h_mass_Glo_pass","",MassBinSta_, MassMinSta_, MassMaxSta_)
        hGlofail = ROOT.TH1D("h_mass_Glo_fail","",MassBinSta_, MassMinSta_, MassMaxSta_)

        hStapass = ROOT.TH1D("h_mass_Sta_pass","",MassBin_, MassMin_, MassMax_)
        hStafail = ROOT.TH1D("h_mass_Sta_fail","",MassBin_, MassMin_, MassMax_)

    elif factorization_scheme == 2:
        hTrkpass = ROOT.TH1D("h_mass_Trk_pass","",MassBinSta_, MassMinSta_, MassMaxSta_)
        hTrkfail = ROOT.TH1D("h_mass_Trk_fail","",MassBinSta_, MassMinSta_, MassMaxSta_)

        # hTrkpass2 = ROOT.TH1D("h_mass_Trk_pass2","",MassBinSta_, MassMinSta_, MassMaxSta_)
    
    byLS_data = byLS_data.loc[(byLS_data['run'] >= int(args.beginRun)) & (byLS_data['run'] < int(args.endRun))]

    recLumi = 0
    firstRun = 0
    lastRun = 0
    df=None
    results = []
    mergeNextRun=False
    log.info(" === Looping over runs... {0} to {1}".format(int(args.beginRun), int(args.endRun)))
    for run, byLS_run in byLS_data.groupby('run', sort=True):
        
        # first and last run of the measurement
        if firstRun == 0:
            firstRun = run
        lastRun = run

        fill = byLS_run.drop_duplicates('fill')['fill'].values[0]
        LSlist = byLS_run['ls'].values.tolist()

        log.info(" === Running Fill {0}".format(fill))
        log.info(" === Running Run {0}".format(run))
        
        eosFile = eosDir+"/"+prefix_dqm+getEra(run)+"*Muon_"+str(run)+"*.root"
        eosFiles = glob.glob(eosFile)
        if len(eosFiles) == 1:
            eosFile = eosFiles[0]
        else:
            log.warning(" === No file or more than one was found! - continue")
            log.warning(" === Was looking for: {}".format(eosFile))            
            continue
        file_ = ROOT.TFile(eosFile,"READ")

        # histograms need to be in same directory so that they can get filled
        hPV.SetDirectory(file_)
        h2HLT.SetDirectory(file_)
        h1HLT.SetDirectory(file_)
        hSelfail.SetDirectory(file_)

        # trees with muon pairs
        tHLT = file_.Get("HLT")
        tSel = file_.Get("Sel")

        if factorization_scheme == 1:
            tGlo = file_.Get("Glo")
            tSta = file_.Get("Sta")

            hSelpass.SetDirectory(file_)

            hGlopass.SetDirectory(file_)
            hGlofail.SetDirectory(file_)

            hStapass.SetDirectory(file_)
            hStafail.SetDirectory(file_)

        elif factorization_scheme == 2:
            tTrk = file_.Get("Trk")

            hTrkpass.SetDirectory(file_)
            hTrkfail.SetDirectory(file_)

            # hTrkpass2.SetDirectory(file_)

        Lumilist = byLS_run.loc[byLS_run['ls'].isin(LSlist)]['recorded(/pb)'].values.tolist()
        ZCountlist = [tHLT.GetEntries("lumiBlock=={0}".format(l)) for l in LSlist]

        log.debug(" === Have lumi secion list {0}".format(LSlist))        
        log.info(" === Looping over measurements...")
        for m, goodLSlist in enumerate(
            get_ls_for_next_measurement(lumisections=LSlist, luminosities=Lumilist, zcounts=ZCountlist, 
                lumiPerMeasurement=LumiPerMeasurement)
        ):
            log.debug(" === Selected lumi section list {0}".format(goodLSlist))

            if measurement is not None and measurement < m:
                break
                    
            # create datafram byLS for measurement
            byLS_m = byLS_run.loc[byLS_run['ls'].isin(goodLSlist)]
            
            ### fill histograms
            file_.cd() # switch to directory where ttrees and histograms are placed

            # define acceptance cuts
            acceptance = " && mass>={0} && mass<{1} && ptTag > {2} && ptProbe > {2} && abs(etaTag) < {3} && abs(etaProbe) < {3}".format(MassMin_, MassMax_, args.ptCut, args.etaCut)
            acceptanceSta = " && mass>={0} && mass<{1} && ptTag > {2} && ptProbe > {2} && abs(etaTag) < {3} && abs(etaProbe) < {3}".format(MassMinSta_, MassMaxSta_, args.ptCut, args.etaCut)

            log.info(" === Fill histograms for measurement {0} ...".format(m))                        
            for iLS in goodLSlist:
                    
                tHLT.Draw("nPV>>+h_PV","lumiBlock=={0}".format(iLS))
                
                n2Before = h2HLT.Integral()
                n1Before = h1HLT.Integral()

                tHLT.Draw("mass>>+h_mass_2HLT_Z",   "pass==2  && lumiBlock=={0} {1}".format(iLS, acceptance))
                tHLT.Draw("mass>>+h_mass_1HLT_Z",   "pass==1  && lumiBlock=={0} {1}".format(iLS, acceptance))
                tSel.Draw("mass>>+h_mass_Sel_fail", "pass==0  && lumiBlock=={0} {1}".format(iLS, acceptance))
                
                if factorization_scheme == 1:
                    tSel.Draw("mass>>+h_mass_Sel_pass", "pass==1  && lumiBlock=={0} {1}".format(iLS, acceptance))

                    tGlo.Draw("mass>>+h_mass_Glo_pass","pass==1 && lumiBlock=={0} {1}".format(iLS, acceptanceSta))                
                    tGlo.Draw("mass>>+h_mass_Glo_fail","pass==0 && lumiBlock=={0} {1}".format(iLS, acceptanceSta))                
 
                    tSta.Draw("mass>>+h_mass_Sta_pass","pass==1 && lumiBlock=={0} {1}".format(iLS, acceptance))   
                    tSta.Draw("mass>>+h_mass_Sta_fail","pass==0 && lumiBlock=={0} {1}".format(iLS, acceptance))  

                elif factorization_scheme == 2:

                    tTrk.Draw("mass>>+h_mass_Trk_pass", "pass>=1  && lumiBlock=={0} {1}".format(iLS, acceptanceSta))   
                    tTrk.Draw("mass>>+h_mass_Trk_fail", "pass==0  && lumiBlock=={0} {1}".format(iLS, acceptanceSta))   

                    # tTrk.Draw("mass>>+h_mass_Trk_pass2", "pass==2  && lumiBlock=={0} {1}".format(iLS, acceptanceSta))   

                n2After = h2HLT.Integral()
                n1After = h1HLT.Integral()
                
                # store the number of 1hlt and 2hlt events in each lumisection
                n2 = n2After - n2Before
                n1 = n1After - n1Before

                byLS_m.loc[byLS_m['ls'] == iLS, 'N2HLT'] = n2
                byLS_m.loc[byLS_m['ls'] == iLS, 'N1HLT'] = n1                    

            if df is None:
                df = byLS_m
            else:
                df = df.append(byLS_m, sort=False)
            
            recLumi = df['recorded(/pb)'].sum()

            log.info(" === Have now recorded lumi = {0}".format(recLumi))            
            log.info(" === Have now {0} | {1} events".format(df['N2HLT'].sum(), h2HLT.Integral()))
            log.info(" === Histograms filled ...")  
            
            # check if upcoming runs make enough data for a measurement
            lumi=0
            mergeNextRun=True
            nextRun = run+1
            for r, dr in byLS_data.groupby('run'):
                if r <= run or r >= int(args.endRun):
                    continue
                # if getFileName(eosDir, r) is None: # check if file of next run exists
                #     continue
                nextRun = r
                LS = dr['ls'].values.tolist()
                lumi += sum(dr.loc[dr['ls'].isin(LS)]['recorded(/pb)'].values)
                if lumi > 0.5 * LumiPerMeasurement:
                    mergeNextRun = False
                    break
            
            mergeNextRun = nextRun < int(args.endRun) and (mergeNextRun or recLumi < 0.5 * LumiPerMeasurement or args.inclusive)            

            if mergeNextRun:
                log.info(" === Merge with next run ... ")
                continue
            
            if firstRun != lastRun:
                outSubDir = outDir + "Run{0}to{1}".format(firstRun,lastRun)
            else:
                outSubDir = outDir + "Run{0}/".format(run)

            log.debug(" === Running measurement {0}".format(m))

            if not args.collect:
                
                if measurement is None or measurement == m:
                    # skip the fit if we look for another measurement
                
                    if not os.path.isdir(outSubDir):
                        os.mkdir(outSubDir)
                    
                    ROOT.set_output(outSubDir)
                    ROOT.set_luminosity(recLumi)
    
                    ROOT.calculateHLTEfficiencyAndYield(h2HLT, h1HLT, m, etaRegionZ, sigModel, bkgModelPass, sigModel, bkgModelPass, hPV, sigTemplates )
                    if factorization_scheme == 1:

                        ROOT.calculateDataEfficiency(hSelpass, hSelfail, m, "Sel",etaRegion, sigModel, bkgModelPass, sigModel, bkgModelFail, 0, sigTemplates )
                        ROOT.set_massRange(MassMinSta_, MassMaxSta_, MassBinSta_)
                        ROOT.calculateDataEfficiency(hGlopass, hGlofail, m, "Glo", etaRegion, sigModel, bkgModelPass, sigModel, bkgModelFail, 0, sigTemplates )
                        ROOT.set_massRange(MassMin_, MassMax_, MassBin_)
                        ROOT.calculateDataEfficiency(hStapass, hStafail, m, "Sta", etaRegion, sigModel, bkgModelPass, sigModel, bkgModelFail, 0, sigTemplates )
                    
                    elif factorization_scheme == 2:
                        ROOT.getZyield(hSelfail, m, "Sel", etaRegion, sigModel, bkgModelFail, 0, sigTemplates, hPV)

                        ROOT.set_massRange(MassMinSta_, MassMaxSta_, MassBinSta_)
                        ROOT.calculateDataEfficiency(hTrkpass, hTrkfail, m, "Trk", etaRegion, sigModel, bkgModelPass, sigModel, bkgModelFail, 0, sigTemplates )
                        
                        # if hTrkpass2.Integral() > 0:
                            # ROOT.getZyield(hTrkpass2, m, "Trk", etaRegion, sigModel, bkgModel, 2, sigTemplates, 0)

                        ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

                    # remove the histogram templates, not needed anymore
                    os.system("rm {0}/histTemplates_*".format(outSubDir))

            cIO = getCorrelation(hPV, mcCorrelations, which="IO")

            cID = getCorrelation(hPV, mcCorrelations, which="ID")

            result = extract_results(outSubDir, m, cIO, cID)
            
            if result:              
                delLumi = df['delivered(/pb)'].sum()
                recLumi = df['recorded(/pb)'].sum()
                # deadtime - fraction during the good lumisections
                deadtime = recLumi / delLumi
                timewindow = len(df) * secPerLS

                # convert time string to datetime format
                beginTime = to_DateTime(df['time'][0], string_format = "mm/dd/yy")
                endTime = to_DateTime(df['time'][-1], string_format = "mm/dd/yy") + datetime.timedelta(seconds=secPerLS)
                
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
                log.info(" === No result - continue")
            
            # prepare for next measurement
            df=None

            # clean the histograms for the next measurement
            hPV.Reset() 
            h2HLT.Reset()
            h1HLT.Reset()
            hSelfail.Reset()

            if factorization_scheme == 1:
                hSelpass.Reset()
                hGlopass.Reset()
                hGlofail.Reset()
                hStapass.Reset()
                hStafail.Reset()
            elif factorization_scheme == 2:
                hTrkpass.Reset()
                hTrkfail.Reset()
                # hTrkpass2.Reset()

        ### prepare for next run
        # keep histograms
        hPV.SetDirectory(0)
        h2HLT.SetDirectory(0)
        h1HLT.SetDirectory(0)
        hSelfail.SetDirectory(0)

        if factorization_scheme == 1:
            hSelpass.SetDirectory(0)
            hGlopass.SetDirectory(0)
            hGlofail.SetDirectory(0)
            hStapass.SetDirectory(0)
            hStafail.SetDirectory(0)
        elif factorization_scheme == 2:
            hTrkpass.SetDirectory(0)
            hTrkfail.SetDirectory(0)  

            # hTrkpass2.SetDirectory(0)

        file_.Close()

        if mergeNextRun:
            continue
        
        if measurement is None or measurement == m:
            ## Write per measurement csv file - one per run
            log.info(" === Writing per Run CSV file")
            results = pd.concat([pd.DataFrame([result]) for result in results], ignore_index=True, sort=False)

            with open(outCSVDir + '/csvfile{0}.csv'.format(run), 'w') as file:
                results.to_csv(file, index=False)

        firstRun = 0
        results = []

    if args.writeSummaryCSV:
        writeSummaryCSV(outCSVDir, writeByLS=False)

    log.info(" ===Done")
