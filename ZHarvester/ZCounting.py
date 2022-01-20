import logging as log
import ROOT
import pandas as pd
import glob
import os
import numpy as np
import json
import pdb
import uncertainties as unc
import gc

from python.utils import writeSummaryCSV

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime, getMCCorrection, unorm, pquad, plinear

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None

# turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

# <<<------------------------
def do_fits(m,
    h2HLTBB_, h2HLTBE_, h2HLTEE_,
    h1HLTBB_, h1HLTBE_, h1HLTEE_,
    fit_one_bins, fit_two_bins,
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
    hPV=0,
):

    if fit_two_bins:
        # Fit BB, BE and EE separately but fit 2HLT and 1HLT each time together
        ROOT.calculateHLTEfficiencyAndYield(h2HLTBB_, h1HLTBB_,
            m, "BB", sigmod_yield, bkgmod_yield, sigmod_yield, bkgmod_yield,
            # cBB_, 
            hPV, sigtempl_yield)
        ROOT.calculateHLTEfficiencyAndYield(h2HLTBE_, h1HLTBE_,
            m, "BE", sigmod_yield, bkgmod_yield, sigmod_yield, bkgmod_yield,
            # cBE_, 
            hPV, sigtempl_yield)
        ROOT.calculateHLTEfficiencyAndYield(h2HLTEE_, h1HLTEE_,
            m, "EE", sigmod_yield, bkgmod_yield, sigmod_yield, bkgmod_yield,
            # cEE_, 
            hPV, sigtempl_yield)

        # fit muon efficiency in two eta bins, barrel (B) and endcap (E)
        ROOT.calculateGloEfficiency(h1TrkBPass, h1TrkBFail, h1StaBFail,
            m, "B", sigmod_trk[0], bkgmod_trk[0], sigmod_trk[1], bkgmod_trk[1],
            0, sigtempl_trk, bkgshape_trk)
        ROOT.calculateGloEfficiency(h1TrkEPass, h1TrkEFail, h1StaEFail,
            m, "E", sigmod_trk[2], bkgmod_trk[2], sigmod_trk[3], bkgmod_trk[3],
            0, sigtempl_trk, bkgshape_trk)

        ROOT.calculateDataEfficiency(h1SITBPass, h1SITBFail,
            m, "Sel", "B", sigmod_sel[0], bkgmod_sel[0], sigmod_sel[1], bkgmod_sel[1],
            0, sigtempl_sel, bkgshape_sel)
        ROOT.calculateDataEfficiency(h1SITEPass, h1SITEFail,
            m, "Sel", "E", sigmod_sel[2], bkgmod_sel[2], sigmod_sel[3], bkgmod_sel[3],
            0, sigtempl_sel, bkgshape_sel)

    if fit_one_bins:
        # fit muon efficiency in one eta bin, inclusive (I)
        # Merge the histograms together
        h1TrkPass = h1TrkBPass.Clone()
        h1TrkFail = h1TrkBFail.Clone()
        h1StaPass = h1StaBPass.Clone()
        h1StaFail = h1StaBFail.Clone()
        h1SITPass = h1SITBPass.Clone()
        h1SITFail = h1SITBFail.Clone()
        
        h1TrkPass.Add(h1TrkEPass)
        h1TrkFail.Add(h1TrkEFail)
        h1StaPass.Add(h1StaEPass)
        h1StaFail.Add(h1StaEFail)
        h1SITPass.Add(h1SITEPass)
        h1SITFail.Add(h1SITEFail)
                        
        h2HLT = h2HLTBB_.Clone()
        h2HLT.Add(h2HLTBE_)
        h2HLT.Add(h2HLTEE_)

        h1HLT = h1HLTBB_.Clone()
        h1HLT.Add(h1HLTBE_)
        h1HLT.Add(h1HLTEE_)       
        
        # --- combined fits to directly extract efficiencies
        # Fit BB, BE and EE separately but fit 1HLT and 2HLT each time together        

        ROOT.calculateHLTEfficiencyAndYield(h2HLT, h1HLT,
            m, "I", sigmod_yield, bkgmod_yield, sigmod_yield, bkgmod_yield,
            # cI_, 
            hPV, sigtempl_yield)        

        ROOT.calculateDataEfficiency(h1SITPass, h1SITFail,
            m, "Sel", "I", sigmod_sel[0], bkgmod_sel[0], sigmod_sel[1], bkgmod_sel[1],
            0, sigtempl_sel, bkgshape_sel)

        ROOT.calculateGloEfficiency(h1TrkPass, h1TrkFail, h1StaFail,
            m, "I", sigmod_trk[0], bkgmod_trk[0], sigmod_trk[1], bkgmod_trk[1],
            0, sigtempl_trk, bkgshape_trk)

def extract_results(directory, m):#, cBB=1, cBE=1, cEE=1, cI=1):
    log.info(" ===Extracting fit results in {0} for {1}".format(directory,m))
    
    res = {}
    
    for eta in ("B","E","I"):
        
        # --- extract results from combined fits
        file_glo = directory+"/workspace_Glo_{0}_{1}.root".format(eta, m)
        if not os.path.isfile(file_glo):
            continue
            
        fGlo = ROOT.TFile(file_glo,"READ")
        rGlo = fGlo.Get("fitResult")
        wGlo = fGlo.Get("workspace")
        effTrk = wGlo.var("effTrk")
        effSta = wGlo.var("effSta")
        (Trkeff, Staeff) = unc.correlated_values_norm(
            [(effTrk.getVal(), effTrk.getError()),(effSta.getVal(), effSta.getError())],
            [[1.,rGlo.correlation("effTrk","effSta")],[rGlo.correlation("effTrk","effSta"),1.]])
        Gloeff = Trkeff*Staeff
        
        fSel = ROOT.TFile(directory+"/workspace_Sel_{0}_{1}.root".format(eta, m),"READ")
        wSel = fSel.Get("workspace")
        Seleff = unc.ufloat(wSel.var("eff").getVal(), wSel.var("eff").getError())
        
        res["Seleff"+eta] = Seleff
        res["Gloeff"+eta] = Gloeff
        res["Trkeff"+eta] = Trkeff
        res["Staeff"+eta] = Staeff
        res["Seleff{0}_chi2".format(eta)] = wSel.arg("chi2").getVal()
        res["Gloeff{0}_chi2".format(eta)] = wGlo.arg("chi2").getVal()
        res["Seleff{0}_chi2pass".format(eta)] = wSel.arg("chi2pass").getVal()
        res["Seleff{0}_chi2fail".format(eta)] = wSel.arg("chi2fail").getVal()
        res["Gloeff{0}_chi2pass".format(eta)] = wGlo.arg("chi2pass").getVal()
        res["Trkeff{0}_chi2fail".format(eta)] = wGlo.arg("chi2failTrk").getVal()
        res["Staeff{0}_chi2fail".format(eta)] = wGlo.arg("chi2failSta").getVal()

        fGlo.Close()
        fGlo.Delete()
        wGlo.Delete()
        rGlo.Delete()
        fSel.Close()
        fSel.Delete()
        wSel.Delete()

    for eta in ("BB","BE","EE","I"):
        
        file_yield = directory+"/workspace_{0}_{1}.root".format(eta, m)
        
        if not os.path.isfile(file_yield):
            continue
        
        f = ROOT.TFile(file_yield,"READ")
        w = f.Get("workspace")
        HLTeff  = unc.ufloat(w.var("eff").getVal(), w.var("eff").getError())
        zRec = unc.ufloat(w.var("Nsig").getVal(), w.var("Nsig").getError())
        NbkgHLT1 = unc.ufloat(w.var("NbkgFail").getVal(), w.var("NbkgFail").getError())
        NbkgHLT2 = unc.ufloat(w.var("NbkgPass").getVal(), w.var("NbkgPass").getError())
        corr  = w.arg("c").getVal()
    
        if eta == "I":
            Zeff = res["SeleffI"] * res["SeleffI"] * res["GloeffI"] * res["GloeffI"]
        else:
            Zeff = res["Seleff"+eta[0]] * res["Seleff"+eta[1]] * res["Gloeff"+eta[0]] * res["Gloeff"+eta[1]]
            
        res["HLTeff"+eta] = HLTeff
        res["c{0}".format(eta)] = corr
        res["zRec"+eta] = zRec
        res["Z{0}eff".format(eta)] = Zeff
        res["zDel"+eta] = zRec / Zeff
        res["NbkgHLT1"+eta] = NbkgHLT1
        res["NbkgHLT2"+eta] = NbkgHLT2
        res["yield{0}_chi2".format(eta)] = w.arg("chi2").getVal()
        res["yield{0}_chi2pass".format(eta)] = w.arg("chi2pass").getVal()
        res["yield{0}_chi2fail".format(eta)] = w.arg("chi2fail").getVal()
        
        f.Close()
        f.Delete()
        w.Delete()

    return res

# <<<------------------------
def add_ZCandidates_perLS(df_byLS, eosFilt):
    # add information to the dataframe about how many events were found in the 2 HLT and 1 HLT category

    file_ = ROOT.TFile(eosFile)
    h2HLTBB_ls = file_.Get("h_mass_2HLTBB_Z")
    h1HLTBB_ls = file_.Get("h_mass_1HLTBB_Z")
    h2HLTBE_ls = file_.Get("h_mass_2HLTBE_Z")
    h1HLTBE_ls = file_.Get("h_mass_1HLTBE_Z")
    h2HLTEE_ls = file_.Get("h_mass_2HLTEE_Z")
    h1HLTEE_ls = file_.Get("h_mass_1HLTEE_Z")

    binLo = h2HLTBB_ls.GetYaxis().FindBin(MassMin_)
    binHi = h2HLTBB_ls.GetYaxis().FindBin(MassMax_)

    for eta, h1, h2 in (
        ("BB", h1HLTBB_ls, h2HLTBB_ls),
        ("BE", h1HLTBE_ls, h2HLTBE_ls),
        ("EE", h1HLTEE_ls, h2HLTEE_ls)):

        df_byLS['N1HLT{}'.format(eta)] = df_byLS['ls'].apply(lambda ls, h=h1: h.ProjectionY("", ls, ls, "e").Integral(binLo, binHi))
        df_byLS['N2HLT{}'.format(eta)] = df_byLS['ls'].apply(lambda ls, h=h2: h.ProjectionY("", ls, ls, "e").Integral(binLo, binHi))

    df_byLS['N1HLTI'] = df_byLS['N1HLTBB'] + df_byLS['N1HLTBE'] + df_byLS['N1HLTEE']
    df_byLS['N2HLTI'] = df_byLS['N2HLTBB'] + df_byLS['N2HLTBE'] + df_byLS['N2HLTEE']

# <<<------------------------
# compute how many Z are delivered for each lumisection
def ls_corrections(
    df_byLS,                # Dataframe with one entry per lumisection
    result_m,               # Result by measurement
    m,
    sort_out_ls = False,
    currentYear = 2017
):    
    if sort_out_ls:
        # sort out lumisections without any Z candidate (maybe trigger was off)
        # consider lumisections where we would expect to have at least any count (> 0.01 /pb we expect 0.01*500 = 5 Z bosons)
        zSum = df_byLS['N1HLTI'] + df_byLS['N2HLTI']
        zeros = (zSum == 0) & (df_byLS['recorded(/pb)'] > 0.01)
        if sum(zeros) > 0:
            log.info("Sort out {0} lumisections in run {1}: {2}".format(sum(zeros), run, df_byLS[zeros]['ls'].values))
            df_byLS = df_byLS[~zeros]

    # --- per measurement result
    df_byLS['time'] = df_byLS['time'].apply(lambda x: to_RootTime(x,currentYear))
    # df_byLS['run'] = np.ones(len(df_byLS),dtype='int') * run
    # df_byLS['fill'] = np.ones(len(df_byLS),dtype='int') * fill

    delLumi_m = df_byLS['delivered(/pb)'].sum()
    recLumi_m = df_byLS['recorded(/pb)'].sum()
    deadtime_m = recLumi_m / delLumi_m
    timeWindow_m = len(df_byLS) * secPerLS

    res = {
        "fill": fill,
        "run": run,
        "measurement": m,
        "tdate_begin": min(df_byLS['time']),
        "tdate_end": max(df_byLS['time']),
        "lumiDel": delLumi_m,
        "lumiRec": recLumi_m,
        "timewindow": timeWindow_m,
        "deadtime": deadtime_m,
        "pileUp": df_byLS['avgpu'].mean(),
    }
        
    # --- per ls dataframe, one per measurement
    for eta in ("BB", "BE", "EE", "I" ):   
        if "zRec"+eta not in result_m.keys():
            continue
                        
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

        df_byLS['Nsig1HLT{}'.format(eta)] = purity1_m.n * df_byLS['N1HLT{}'.format(eta)]
        df_byLS['Nsig2HLT{}'.format(eta)] = purity2_m.n * df_byLS['N2HLT{}'.format(eta)]

        # eff_HLT = 2*n2 /(c*(n1 + 2*n2))
        df_byLS['HLTeff{}_mc'.format(eta)] = 2*df_byLS['Nsig2HLT{}'.format(eta)] / (corr_m * (df_byLS['Nsig1HLT{}'.format(eta)] + 2*df_byLS['Nsig2HLT{}'.format(eta)]))

        # nZ = (n2 + n1/2.)**2 * c / n2 ( multiplication with c comes later)
        df_byLS['zRec{}_mc'.format(eta)] = (df_byLS['Nsig2HLT{}'.format(eta)] + df_byLS['Nsig1HLT{}'.format(eta)]/2.)**2 / df_byLS['Nsig2HLT{}'.format(eta)]
        
        # set all infinite values to 0
        df_byLS.loc[~np.isfinite(df_byLS['zRec{}_mc'.format(eta)]), 'zRec{}_mc'.format(eta)] = 0

        df_byLS['zRec{}_mcUp'.format(eta)]   = df_byLS['zRec{}_mc'.format(eta)] * (1.3 * (corr_m - 1) + 1)
        df_byLS['zRec{}_mcDown'.format(eta)] = df_byLS['zRec{}_mc'.format(eta)] * (0.7 * (corr_m - 1) + 1)
        df_byLS['zRec{}_mc'.format(eta)] = df_byLS['zRec{}_mc'.format(eta)] * corr_m

        df_byLS['zDel{}_mc'.format(eta)] = df_byLS['zRec{}_mc'.format(eta)] / zEff_m.n
        df_byLS['zDel{}_mcUp'.format(eta)] = df_byLS['zRec{}_mcUp'.format(eta)] / zEff_m.n
        df_byLS['zDel{}_mcDown'.format(eta)] = df_byLS['zRec{}_mcDown'.format(eta)] / zEff_m.n

        res["zRec{0}_mc".format(eta)] = df_byLS['zRec{0}_mc'.format(eta)].sum()
        res["zDel{0}_mc".format(eta)] = df_byLS['zDel{0}_mc'.format(eta)].sum()
        res["zDel{0}_mcUp".format(eta)] = df_byLS['zDel{0}_mcUp'.format(eta)].sum()
        res["zDel{0}_mcDown".format(eta)] = df_byLS['zDel{0}_mcDown'.format(eta)].sum()
        # res["Z{0}eff".format(eta)] = df_byLS['eff{0}'.format(eta)].values.mean(),
        
        res['zDel{}_Up'.format(eta)]   = (NsigHLT2_m + NsigHLT1_m / 2.)**2 * (1.3 * (corr_m - 1) + 1) / NsigHLT2_m / zEff_m
        res['zDel{}_Down'.format(eta)] = (NsigHLT2_m + NsigHLT1_m / 2.)**2 * (0.7 * (corr_m - 1) + 1) / NsigHLT2_m / zEff_m  
        res['zDel{}_Up'.format(eta)]   = res['zDel{}_Up'.format(eta)].n
        res['zDel{}_Down'.format(eta)] = res['zDel{}_Down'.format(eta)].n
        res["zRec{0}".format(eta)] = result_m['zRec{0}'.format(eta)].n
        res["zDel{0}".format(eta)] = result_m['zDel{0}'.format(eta)].n
        

    with open(outCSVDir + 'csvfile{0}_{1}.csv'.format(run, m), 'w') as file:
        df_byLS.to_csv(file, index=False)
    
    return res


# ------------------------------------------------------------------------------
def load_histo(name_, fileName_, lumisections_=[0,], prefix_="", run_=0, suffix_="", pileup=False):
    file_ = ROOT.TFile(fileName_)
    log.debug("===Load histogram {0}".format(name_))
    
    h_X_ls = file_.Get("{0}{1}".format(prefix_,name_)).Clone("{0}{1}{2}".format(prefix_,name_, suffix_))
    h_X_ls.SetDirectory(0)
    h_X = h_X_ls.ProjectionY("h_tmp_{0}_{1}_{2}".format(name_, run_, suffix_), lumisections_[0], lumisections_[0], "e")
    
    for ls in lumisections_[1:]:
        h_X.Add(h_X_ls.ProjectionY("h_tmp_{0}_{1}_{2}_{3}".format(name_, run_, ls, suffix_), ls, ls, "e"))

    log.debug("===Integral = {0}".format(h_X.Integral()))
    
    if pileup:
        h_X.SetDirectory(0)
        return h_X
        
    # create new histogram in correct bin range
    hNew = ROOT.TH1D("h_mass_{0}_{1}_{2}_{3}".format(
        name_, run_, lumisections_[0], suffix_), "",MassBin_, MassMin_, MassMax_)
    log.debug("===Create new histogram in range m = [{0},{1}]".format(MassMin_, MassMax_))
    
    for ibin in range(0, h_X.GetNbinsX()+1):
        binCenter = h_X.GetBinCenter(ibin)
        if binCenter < MassMin_:
            continue
        elif binCenter > MassMax_:
            break
        else:
            content = h_X.GetBinContent(ibin)
            newBin = hNew.FindBin(binCenter)
            hNew.SetBinContent(newBin, hNew.GetBinContent(newBin) + content)
    
    log.debug("===Integral new = {0}".format(hNew.Integral()))
    
    hNew.SetDirectory(0)

    return hNew


# ------------------------------------------------------------------------------
## Get root file with histograms
def getFileName(directory, run):
    # check if run was processed already
    log.info(" ===Checking input DQMIO.root file...")
    eosFileList = glob.glob(directory + '/*/*' + str(run) + '*root')
    if not len(eosFileList) > 0:
        log.info(" The file does not yet exist for run: " + str(run))
        return None
    elif len(eosFileList) > 1:
        log.info(" Multiple files found for run: " + str(run))
        return None
    else:
        return eosFileList[0]

# ------------------------------------------------------------------------------
## load input by lumisection CSV file, convert it and return the data 
def load_input_csv(byLS_data):
    log.info(" Loading input byls csv...")
    byLS_file = open(str(byLS_filename))
    byLS_lines = byLS_file.readlines()
    byLS_data = pd.read_csv(byLS_filename, sep=',', low_memory=False,
        skiprows=lambda x: byLS_lines[x].startswith('#') and not byLS_lines[x].startswith('#run'))
        
    log.info(" formatting csv file...")    # formatting the csv
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

    # if there are multiple entries of the same ls (for example from different triggers), 
    #   only keep the one with the highest luminosity.
    byLS_data = byLS_data.sort_values(['fill', 'run', 'ls', 'delivered(/pb)', 'recorded(/pb)'])
    byLS_data = byLS_data.drop_duplicates(['fill', 'run', 'ls'])
    
    return byLS_data

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
                        help="Directory to root file for signal template ('None': take analytic function) ")
    parser.add_argument("--bkgTemplates", default="None", type=str,
                        help="Directory to root file for bkg template ('None': take analytic function, 'altBkg': take alternative analytic function) ")
    parser.add_argument('--ptCut', type=float, default=30.,
                        help='specify lower pt cut on tag and probe muons')
    parser.add_argument('--etaCut', type=float, default=2.4,
                        help='specify upper |eta| cut on tag and probe muons')
    parser.add_argument('--mass', nargs=3, metavar=('LOW', 'HIGH', 'NUMBER'), default=(70,250,360), type=int,
                        help='specify mass range for tag and probe muon pairs')
    parser.add_argument('--LumiPerMeasurement', default=20, type=float,
                        help='specify amount of luminosity per measurement in pb-1')
    parser.add_argument('--inclusive', default=False, action="store_true",
                        help='specify whether or not to do an inclusive fit of the specified runs')
    parser.add_argument('--collect', default=False, action="store_true",
                        help='specify whether or not to run the fits or just collect the results')
    parser.add_argument("-o", "--dirOut", help="where to store the output files", default="./")

    # Fit the efficiency in two bins: barrel (B) and endcap (E)
    fit_two_bins = True #True
    # Fit the efficiency in one bin: inclusive (I)
    fit_one_bins = True

    args = parser.parse_args()
    if args.verbose:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
    else:
        log.basicConfig(format="%(levelname)s: %(message)s", level=log.INFO)

    ########################################
    # link to resouces
    prefix_dqm="" #"DQMData/Run {0}/ZCounting/Run summary/Histograms/".format(run)
    resPath = cmsswbase + "/src/ZCounting/ZHarvester/res/"
    if( args.beginRun >= 272007 and args.beginRun < 278808
        # there is an overlap for 2016 F in runs with pre and post VFP settings
        and args.beginRun not in [278769, 278801, 278802, 278803, 278804, 278805, 278808]
    ): # 2016 pre VFP
        byLsCSV          = resPath+"/FillByLs_2016.csv"
        sigTemplates     = resPath+"/sigTemplates/V14_IsoMu24/ZSignalTemplate-V14-Summer16preVFP-DYJetsToLL_M_50_LO.root"
        bkgTemplates     = resPath+"/qcdTemplates/Run2016preVFP.root"
        era = "2016preVFP"
    elif args.beginRun < 294645:    # 2016 post VFP
        byLsCSV          = resPath+"/FillByLs_2016.csv"
        sigTemplates     = resPath+"/sigTemplates/V14_IsoMu24/ZSignalTemplate-V14-Summer16postVFP-DYJetsToLL_M_50_LO.root"
        bkgTemplates     = resPath+"/qcdTemplates/Run2016postVFP.root"
        era = "2016postVFP"
    elif args.beginRun > 297020 and args.beginRun < 306828:     # 2017
        byLsCSV          = resPath+"/FillByLs_2017.csv"
        sigTemplates     = resPath+"/sigTemplates/V14_IsoMu27/ZSignalTemplate-V14-Fall17-DYJetsToLL_M_50_LO.root"
        bkgTemplates     = resPath+"/qcdTemplates/Run2017.root"
        era = "2017"
    elif args.beginRun >= 306926 and args.beginRun < 307083:    # 2017 H
        byLsCSV          = resPath+"/FillByLs_2017_lowPU.csv"
        sigTemplates     = resPath+"/sigTemplates/V14_IsoMu24/ZSignalTemplate-V14-Fall17-DYJetsToLL_M_50_LO.root"
        era = "2017H"
    elif args.beginRun >= 315252:                               # 2018
        byLsCSV          = resPath+"/FillByLs_2018.csv"
        sigTemplates     = resPath+"/sigTemplates/V14_IsoMu24/ZSignalTemplate-V14-Autumn18-DYJetsToLL_M_50_LO.root"
        bkgTemplates     = resPath+"/qcdTemplates/Run2018.root"
        era = "2018"
    else:
        byLsCSV             = None
        sigTemplates        = None

    eosDir           = args.dirDQM
    byLsCSV          = byLsCSV          if args.byLsCSV       == "default"   else args.byLsCSV
    bkgTemplates     = bkgTemplates     if args.bkgTemplates  == "default"   else args.bkgTemplates
    measurement       = args.measurement

    log.info("----------------------------------")
    log.info("Use eosDir:              {0}".format(eosDir))
    log.info("Use byLsCSV:             {0}".format(byLsCSV))
    log.info("Use sigTemplates:        {0}".format(sigTemplates))
    log.info("Use bkgTemplates:        {0}".format(bkgTemplates))
    # log.info("Use correlationsFile:    {0}".format(correlationsFile))
    log.info("Mass range from:         {0} to {1}".format(*args.mass))
    log.info("Lumi per Measurement:    {0}".format(args.LumiPerMeasurement))
    log.info("----------------------------------")

    ########## Input configuration ##########
    # ByLS csv inputs generated by testBRIL.sh
    byLS_filelist = glob.glob(byLsCSV)
    byLS_filelist.sort(key=os.path.getmtime)
    byLS_filename = byLS_filelist[-1]
    log.info(" The brilcalc csv file: " + str(byLS_filename))

    # inputs: to build MC*Gaussian template for efficiency fitting
    fit_options = {
        # we always need signal template for calculation of correlation coefficient (c)
        'sigtempl_yield':sigTemplates,
        'sigtempl_sel':sigTemplates,
        'sigtempl_trk':sigTemplates,
        'sigtempl_sta':sigTemplates,
    }
    
    if args.sigTemplates == "default":
        _optns = {
            'sigmod_yield':2,
            'sigmod_sel':[2,2,2,2],
            'sigmod_trk':[2,2,2,2],
            'sigmod_sta':[2,2,2,2],
            }
        fit_options.update(_optns)

    if bkgTemplates == "altBkg":
        _optns = {
            "bkgmod_yield": 4,
            "bkgmod_sel": [4,4,4,4],
            "bkgmod_trk": [4,4,4,4],
            "bkgmod_sta": [4,4,4,4],
            }
        fit_options.update(_optns)
    elif bkgTemplates != "None":
        _optns = {
            "bkgmod_yield": 7,
            "bkgshape_yield": bkgTemplates,
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
            log.warning(": directory already exists ...")

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
        log.info(" Loading C marco...")
        # load functions for fitting
        ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(
            __file__)) + "/calculateDataEfficiency.C")

        ROOT.set_massRange(MassMin_, MassMax_, MassBin_)

        ROOT.set_ptCut(args.ptCut)
        ROOT.set_etaCut(args.etaCut)
    
    byLS_data = load_input_csv(byLS_filename)

    #####################################

    # For fitting
    hPV = ROOT.TH1D("h_PV","", 75, -0.5, 74.5)
    
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
    
    recLumi = 0
    firstRun = 0
    lastRun = 0
    df=None
    results = []
    mergeNextRun=False
    log.info(" ===Looping over runs... {0} to {1}".format(int(args.beginRun), int(args.endRun)))
    for run, data_run in byLS_data.groupby('run'):

        if run < int(args.beginRun) or run >= int(args.endRun):
            continue
        
        # first and last run of the measurement
        if firstRun == 0:
            firstRun = run
        lastRun = run

        fill = data_run.drop_duplicates('fill')['fill'].values[0]
        LSlist = data_run.query('ls <= {0}'.format(maximumLS))['ls'].values.tolist()

        # Consider only runs with at least 5 LS
        if len(LSlist) <=5:
            log.info(" ===Skip Run {0} since it only has {1} LS".format(run, len(LSlist)))
            continue

        log.info(" ===Running Fill {0}".format(fill))
        log.info(" ===Running Run {0}".format(run))

        eosFile = getFileName(eosDir, run)      # >>> histograms binned in mass
        if eosFile is None:
            log.info("Histogram file not found, continue")
            continue

        log.info(" ===Looping over measurements...")
        m = 0
        while len(LSlist) > 0:  # begin next measurement "m"
            log.debug(" ===Have lumi secion list {0}".format(LSlist))

            # merge data to one measuement if remaining luminosity is too less for two measuements
            mergeMeasurements = sum(
                data_run.loc[data_run['ls'].isin(LSlist)]['recorded(/pb)'].values) < 1.5 * LumiPerMeasurement

            # produce goodLSlist with ls that are used for one measurement
            goodLSlist = []
            while len(LSlist) > 0:
                goodLSlist.append(LSlist[0])
                recLumi += (data_run[data_run['ls'] == LSlist[0]]['recorded(/pb)'].values)[0]
                del LSlist[0]
                # if we have enough collected lumisections
                if not mergeMeasurements and recLumi >= LumiPerMeasurement:# or len(goodLSlist) >= LSperMeasurement):
                    break
            
            # create dataframe with one entry per LS
            new_df = data_run.loc[data_run['ls'].isin(goodLSlist)]
            # add info about Z candidates per LS to byLs dataframe
            add_ZCandidates_perLS(new_df, eosFile)                    
            if df is None:
                df = new_df
            else:
                df = df.append(new_df, sort=False)

            log.debug(" ===Selected lumi secion list {0}".format(goodLSlist))

            ### load histograms

            def load(name_):
                return load_histo(name_, eosFile, goodLSlist, prefix_dqm, run, "new")
            
            hPV.Add(load_histo("h_nPV", eosFile, goodLSlist, prefix_dqm, run, "new", pileup=True))
            
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
            
            # check if upcoming runs make enough data for a measurement
            lumi=0
            mergeNextRun=True
            nextRun = run+1
            for r, dr in byLS_data.groupby('run'):
                if r <= run or r >= int(args.endRun):
                    continue
                if getFileName(eosDir, r) is None: # check if file of next run exists
                    continue
                nextRun = r
                LS = dr.query('ls <= {0}'.format(maximumLS))['ls'].values.tolist()
                lumi += sum(dr.loc[dr['ls'].isin(LS)]['recorded(/pb)'].values)
                if lumi > 0.5 * LumiPerMeasurement:
                    mergeNextRun = False
                    break
            
            mergeNextRun = nextRun < int(args.endRun) and (mergeNextRun or recLumi < 0.5 * LumiPerMeasurement or args.inclusive)            
            
            log.debug(" === Have now recorded lumi = {0}".format(recLumi))
            if mergeNextRun:
                continue
            
            if firstRun != lastRun:
                outSubDir = outDir + "Run{0}to{1}".format(firstRun,lastRun)
            else:
                outSubDir = outDir + "Run{0}/".format(run)

            log.debug(" ===Running measurement {0}".format(m))
            
            if not args.collect:
                
                if measurement is None or measurement == m:
                    # skip the fit if we look for another measurement
                
                    if not os.path.isdir(outSubDir):
                        os.mkdir(outSubDir)
                    
                    ROOT.set_output(outSubDir)
                    ROOT.set_luminosity(recLumi)
                        
                    do_fits(m,
                        h2HLTBB, h2HLTBE, h2HLTEE,
                        h1HLTBB, h1HLTBE, h1HLTEE,
                        fit_one_bins, fit_two_bins,
                        hSITpassB, hSITfailB, hSITpassE, hSITfailE,
                        hTrkpassB, hTrkfailB, hTrkpassE, hTrkfailE,
                        hStapassB, hStafailB, hStapassE, hStafailE,
                        hPV=hPV,
                        **fit_options
                        )
                # remove the PV dependent templates 
                os.system('rm {0}/histTemplates*'.format(outSubDir))
                
                # clean the histograms for the next measurement
                h2HLTBB.Reset()
                h2HLTBE.Reset()
                h2HLTEE.Reset()
                h1HLTBB.Reset()
                h1HLTBE.Reset()
                h1HLTEE.Reset()

                hSITpassB.Reset()
                hSITfailB.Reset()
                hSITpassE.Reset()
                hSITfailE.Reset()
                hTrkpassB.Reset()
                hTrkfailB.Reset()
                hTrkpassE.Reset()
                hTrkfailE.Reset()
                hStapassB.Reset()
                hStafailB.Reset()
                hStapassE.Reset()
                hStafailE.Reset()
            
            result = extract_results(outSubDir, m)#, cBB, cBE, cEE, cI)

            if result:
                result.update(ls_corrections(df, result, m,
                    sort_out_ls = (era!="2017H")    #do not sort out any lumisection if we fit low pileup
                    ))
                results.append(result)
            else:
                log.info(" === No result - continue")
            
            # prepare for next measurement
            recLumi = 0     # reset lumi count
            df=None
            m += 1

        if mergeNextRun:
            continue

        ## Write per measurement csv file - one per run
        log.info(" ===Writing per Run CSV file")
        results = pd.concat([pd.DataFrame([result]) for result in results], ignore_index=True, sort=False)

        with open(outCSVDir + '/csvfile{0}.csv'.format(run), 'w') as file:
            results.to_csv(file, index=False)

        results = []


    if args.writeSummaryCSV:
        writeSummaryCSV(outCSVDir)

    log.info(" ===Done")
