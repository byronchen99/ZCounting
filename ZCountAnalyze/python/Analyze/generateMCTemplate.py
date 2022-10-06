#!/usr/bin/env python3
#
################################################################################
#
# create root ttrees from DY MC used as templates in fits
#
################################################################################

import uproot 
import awkward as ak
import numpy as np
import pdb
import vector
import hist
from utils import get_masses
import os

import argparse
parser = argparse.ArgumentParser(prog='./histogramize')
parser.add_argument(
    '-i', '--input', nargs='+',
    help='specify input ntuple root files'
)
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output directory'
)
args = parser.parse_args()


### acceptance cuts
ptCut = 27
etaCut = 2.4
mass_lo = 66
mass_hi = 116

muonMass = 0.105658369

### resources
treename = "zcounting/tree"

filenames = [f+":"+treename for f in args.input]
dirOut = args.output

if not os.path.isdir(dirOut):
    print(f"create output directory {dirOut}")
    os.mkdir(dirOut)


branches_muon = ["Muon_pt", "Muon_eta", "Muon_phi", "Muon_ID", "Muon_charge", "Muon_triggerBits"]
branches_event = ["nPV", "nPU", "eventweight", "decayMode",
    "lepton_genPt", "lepton_genEta", "lepton_genPhi",
    "antiLepton_genPt", "antiLepton_genEta", "antiLepton_genPhi"
    ]
aliases = {

}

hists = dict()
# process the events in batches to avoid using too much memory
for batch in uproot.iterate(filenames, step_size="100 MB", library="ak", 
    filter_name=branches_muon+branches_event, aliases=aliases
):
    if len(batch) == 0:
        continue
    print(f"Have a batch with {len(batch)} events")    
    
    muons = batch[branches_muon]
    events = batch[branches_event]
    
    # select muons within the acceptance
    muons = muons[(muons["Muon_pt"] > ptCut) & (abs(muons["Muon_eta"]) < etaCut)]

    # select events with at least two muons from which at least one passes ID and is an HLT muon 
    mask = (ak.num(muons["Muon_pt"]) >= 2) 
    mask = mask & (ak.num(muons[(muons["Muon_ID"]>=4) & (muons["Muon_triggerBits"]&1 != 0)]["Muon_pt"]) >= 1)

    # set muon rest mass
    muons["Muon_mass"] = muons["Muon_pt"]*0+muonMass
    
    # make pairs of two muons
    pairs = ak.combinations(muons[["Muon_charge", "Muon_ID", "Muon_triggerBits", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"]], 2)
    
    # require pairs with opposite charge
    lefts, rights = ak.unzip(pairs["Muon_charge"])
    pairs = pairs[lefts*rights == -1]
    
    # of each pair, require at least one muon pass the ID and HLT_IsoMu24 - the tag muon
    l_HLT, r_HLT = ak.unzip(pairs["Muon_triggerBits"])
    l_ID, r_ID = ak.unzip(pairs["Muon_ID"])
    pairs = pairs[((l_HLT&1 != 0) & (l_ID>=4)) | ((r_HLT&1 != 0) & (r_ID>=4))]    
    
    # for further selections
    l_HLT, r_HLT = ak.unzip(pairs["Muon_triggerBits"])
    l_ID, r_ID = ak.unzip(pairs["Muon_ID"])

    p_hlt1 = pairs[(l_ID >= 4) & (r_ID >= 4) & ((l_HLT&1) != (r_HLT&1))] # require both muons to pass the ID and one muon pass the trigger
    p_hlt2 = pairs[(l_ID >= 4) & (r_ID >= 4) & (l_HLT&1 != 0) & (r_HLT&1 != 0)]  # require both muons to pass the ID and both muons pass the trigger
    p_sel_fail = pairs[(l_ID >= 3) & (r_ID >= 3) & ((l_ID >= 4) != (r_ID >= 4))] # require both muons to be isGlobalMuon but one muon fails ID
    
    # loop over pair collections with names to store them in dictionary
    for name, pair in ( 
        ("hlt1", p_hlt1),
        ("hlt2", p_hlt2), 
        ("sel_fail", p_sel_fail) 
    ):
        # get masses for each pair    
        events[f"m_{name}"] = get_masses(pair, mass_lo, mass_hi)
        
        mask = (ak.num(events[f"m_{name}"]) >= 1)

        xx = events[mask]["lumiBlock"]
        yy = events[mask][f"m_{name}"]
        
        # pair multiplicities in each event 
        # counts = ak.num(events[mask][f"m_{n}"])
        
        # bring the event by event array into the same dimension
        if len(xx) > 0:
            xx = ak.unflatten(xx,counts=1)
        
        # bring the arrays into the same form
        xx, yy = ak.broadcast_arrays(xx,yy)
        
        xx = ak.flatten(xx).to_numpy()
        yy = ak.flatten(yy).to_numpy()
        
        hist = np.histogram2d(xx, yy, bins=[lumi_nBins, mass_nBins], range=((lumi_lo,lumi_hi), (mass_lo, mass_hi)))
        
        histname = f"h_mass_{name}"
        if histname in hists.keys():
            hists[histname] = (hist[0] + hists[histname][0], hist[1], hist[2])
        else:
            hists[histname] = hist

# open output root file
print(f"Write out results in `output.root`")
output = uproot.recreate(f"{dirOut}/output.root")

# write the histograms out
for name, hist in hists.items():
    output[name] = hist




