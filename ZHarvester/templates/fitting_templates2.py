import uproot
import numpy as np
import awkward as ak
import vector
from hist import Hist
import coffea
from coffea import lookup_tools

import argparse
import h5py

def _get_rochester_corr(rochester_filename):
    path = rochester_filename
    
    rochester_data = lookup_tools.txt_converters.convert_rochester_file(
        path, loaduncs=True
    )
    rochester = lookup_tools.rochester_lookup.rochester_lookup(
        rochester_data
    )
    return rochester


def generate_templates(data_reco, rochester):

    # get rochester corrected pt
    for inner_type in ["Muon", "Track"]:
        data_reco["correct_charge_genPt"] = (data_reco["{0}_charge".format(inner_type)] + 1)/2 * data_reco["antiMuon_genPt"] - (data_reco["{0}_charge".format(inner_type)] - 1)/2 * data_reco["muon_genPt"]


        scale = rochester.kScaleDT(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)])
        spread = rochester.kSpreadMC(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)], data_reco["correct_charge_genPt"])

        scale_error = rochester.kScaleDTerror(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)])
        spread_error = rochester.kSpreadMCerror(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)], data_reco["correct_charge_genPt"])


        data_reco["{0}_pt_RC_nominal".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * spread
        data_reco["{0}_pt_RC_plus_scale".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * (scale + scale_error) * spread
        data_reco["{0}_pt_RC_minus_scale".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * (scale - scale_error) * spread
        data_reco["{0}_pt_RC_plus_spread".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * (spread + spread_error)
        data_reco["{0}_pt_RC_minus_spread".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * (spread - spread_error)


    template_types = ["HLT2", "HLT1", "HLT0",
     "ID2", "ID1", "ID0", "IDfail"
     "Glo2", "Glo1", "Glo0", "Glofail",
     "Sta2", "Sta1", "Sta0", "Stafail"]

    for template
        