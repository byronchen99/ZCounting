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