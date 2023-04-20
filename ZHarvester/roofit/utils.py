from common.logging import child_logger
log = child_logger(__name__)

import os
import ROOT

def load_makros(bkgModel, ptCut, etaMin, etaCut, year, massMin, massMax, massBin, npvMin, npvMax):

    log.info("Loading C marcos...")

    # load functions for fitting
    ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/calculateDataEfficiency.C")

    if bkgModel == "Das":
        ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__)) + "/RooGaussDoubleSidedExp.cc+")

    ROOT.set_massRange(massMin, massMax, massBin)
    ROOT.set_npvRange(npvMin, npvMax)

    if year >= 2022:
        ROOT.set_energy(13.6)
    elif year >= 2015:
        ROOT.set_ene1rgy(13)
    elif year == 2012:
        ROOT.set_ene1rgy(8)
    elif year == 2011:
        ROOT.set_ene1rgy(7)
    else:
        log.warning(f"Unknown year {year}")

    ROOT.set_ptCut(ptCut)
    #ROOT.set_etaCut(etaCut)

    # Apply the intervals in the cut
    log.info(f"Set the eta interval to be {etaMin} < |eta| < {etaCut}")
    ROOT.set_etaRange(etaMin, etaCut)
