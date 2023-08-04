import uproot
import coffea
from coffea import lookup_tools

import awkward as ak
import numpy as np

import argparse

import pdb


def get_corrections(data_reco, rochester):

    for inner_type in ["Muon", "Track"]:
        data_reco["correct_charge_genPt"] = (data_reco["{0}_charge".format(inner_type)] + 1)/2 * data_reco["antiMuon_genPt"] - (data_reco["{0}_charge".format(inner_type)] - 1)/2 * data_reco["muon_genPt"]
    
        # rochester correction for spread and scale
        scale = rochester.kScaleDT(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)])
        spread = rochester.kSpreadMC(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)], data_reco["correct_charge_genPt"])

        scale_error = rochester.kScaleDTerror(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)])
        spread_error = rochester.kSpreadMCerror(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)], data_reco["correct_charge_genPt"])

        # corrections to pt
        data_reco["{0}_pt_RC_nominal".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * spread
        data_reco["{0}_pt_RC_plus_scale".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * (scale + scale_error) * spread
        data_reco["{0}_pt_RC_minus_scale".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * (scale - scale_error) * spread
        data_reco["{0}_pt_RC_plus_spread".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * (spread + spread_error)
        data_reco["{0}_pt_RC_minus_spread".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * (spread - spread_error)
    
    return data_reco


parser = argparse.ArgumentParser()
parser.add_argument("mc_data_file")
parser.add_argument("rochester_file")
args = parser.parse_args()


filename_ntuples = args.mc_data_file
rochester_filename = args.rochester_file

def _get_rochester_corr(rochester_filename):
    path = rochester_filename
    
    rochester_data = lookup_tools.txt_converters.convert_rochester_file(
        path, loaduncs=True
    )
    rochester = lookup_tools.rochester_lookup.rochester_lookup(
        rochester_data
    )
    return rochester

rochester = _get_rochester_corr(rochester_filename)

# f1 = uproot.open(filename_ntuples)

treename = "zcounting/tree"
# data = uproot.open(filename_ntuples+":"+treename)

columns = ["muon_genPt","muon_genEta", "muon_genPhi",
 "antiMuon_genPt","antiMuon_genEta","antiMuon_genPhi",
 "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge",
 "Muon_ID", "Muon_triggerBits", "nPV",
 "Muon_useUpdated", "Muon_ptStaUpd", "Muon_etaStaUpd", "Muon_phiStaUpd", "Muon_chargeStaUpd",
 "Muon_ptStaReg", "Muon_etaStaReg", "Muon_phiStaReg", "Muon_chargeStaReg",
 "Track_pt", "Track_eta", "Track_phi", "Track_charge",
 "Muon_ptTrk", "Muon_etaTrk", "Muon_phiTrk", "Muon_chargeTrk",
 "Track_nPixelHits", "Track_nTrackerLayers", "Track_trackAlgo",
 'Muon_nPixelHits', 'Muon_nTrackerLayers','Muon_trackAlgo', "decayMode", "nMuon"]


new_columns = ["Muon_pt_RC_nominal", "Muon_pt_RC_plus_scale", "Muon_pt_RC_minus_scale", "Muon_pt_RC_plus_spread", "Muon_pt_RC_minus_spread",
 "Track_pt_RC_nominal", "Track_pt_RC_plus_scale", "Track_pt_RC_minus_scale", "Track_pt_RC_plus_spread", "Track_pt_RC_minus_spread"]

keys = columns + new_columns

new_data_reco = {}
for key in keys:
    new_data_reco[key] = np.array([])

i=0

for array in uproot.iterate(filename_ntuples+":"+treename, expressions = columns, step_size=100000):
    data_reco = array[columns]
    data_reco1 = data_reco[(data_reco["decayMode"]==13)]
    # pdb.set_trace()
    data_reco2 = get_corrections(data_reco1, rochester)
    i=i+1
    print(i)
    
    for key in keys:
        new_data_reco[key] = ak.concatenate((new_data_reco[key], data_reco2[key]))
    
    break


output = uproot.create("MC_with_corrections.root")


for key in keys:
    output[key] = new_data_reco[key]