import awkward as ak
import numpy as np
from hist import Hist
import vector
import uproot
import coffea
from coffea import lookup_tools

import h5py
import argparse
import pdb


def _get_rochester_corr(rochester_filename):
    path = rochester_filename
    
    rochester_data = lookup_tools.txt_converters.convert_rochester_file(
        path, loaduncs=True
    )
    rochester = lookup_tools.rochester_lookup.rochester_lookup(
        rochester_data
    )
    return rochester


# get rochester corrected values of pt
def get_corrections(data_reco, rochester):
    for inner_type in ["Muon", "Track"]:
        data_reco["correct_charge_genPt"] = (data_reco["{0}_charge".format(inner_type)] + 1)/2 * data_reco["antiLepton_genPt"] - (data_reco["{0}_charge".format(inner_type)] - 1)/2 * data_reco["lepton_genPt"]

        scale = rochester.kScaleDT(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)])
        spread = rochester.kSpreadMC(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)], data_reco["correct_charge_genPt"])

        scale_error = rochester.kScaleDTerror(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)])
        spread_error = rochester.kSpreadMCerror(data_reco["{0}_charge".format(inner_type)], data_reco["{0}_pt".format(inner_type)], data_reco["{0}_eta".format(inner_type)], data_reco["{0}_phi".format(inner_type)], data_reco["correct_charge_genPt"])

        data_reco["{0}_pt_RC_none".format(inner_type)] = data_reco["{0}_pt".format(inner_type)]
        data_reco["{0}_pt_RC_nominal".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * spread
        data_reco["{0}_pt_RC_plus_scale".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * (scale + scale_error) * spread
        data_reco["{0}_pt_RC_minus_scale".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * (scale - scale_error) * spread
        data_reco["{0}_pt_RC_plus_spread".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * (spread + spread_error)
        data_reco["{0}_pt_RC_minus_spread".format(inner_type)] = data_reco["{0}_pt".format(inner_type)] * scale * (spread - spread_error)
    
    return data_reco

# define similiar for outer type correction (to be completed)


def get_inner_variables_with_trk(data_reco, correction_type):
    for variable in ["eta", "phi", "charge"]:
        data_reco["inner_{0}".format(variable)] = ak.concatenate((data_reco["Muon_{0}".format(variable)], data_reco["Track_{0}".format(variable)]), axis = 1)
    data_reco["inner_pt_RC_{0}".format(correction_type)] = ak.concatenate((data_reco["Muon_pt_RC_{0}".format(correction_type)], data_reco["Track_pt_RC_{0}".format(correction_type)]), axis = 1)
    data_reco["inner_pt"] = ak.concatenate((data_reco["Muon_pt"], data_reco["Track_pt"]), axis = 1)
    
    for variable in ["nPixelHits", "nTrackerLayers", "trackAlgo"]:
        data_reco["{0}".format(variable)] = ak.concatenate((data_reco["Muon_{0}".format(variable)], data_reco["Track_{0}".format(variable)]), axis = 1)
    
    for variable in ["HLT_trigger_pass", "Muon_ID"]:
        criteria_muons = data_reco["{0}".format(variable)]
        criteria_tracks = ak.ones_like(data_reco["Track_pt"])==0
        data_reco["inner_{0}".format(variable)] = ak.concatenate((criteria_muons, criteria_tracks), axis = 1)
    
    return data_reco

def get_inner_variables_without_trk(data_reco, correction_type):
    for variable in ["eta", "phi", "charge"]:
        data_reco["inner_{0}".format(variable)] = data_reco["Muon_{0}".format(variable)]
    data_reco["inner_pt_RC_{0}".format(correction_type)] = data_reco["Muon_pt_RC_{0}".format(correction_type)]
    data_reco["inner_pt"] = data_reco["Muon_pt"]
    
    for variable in ["HLT_trigger_pass", "Muon_ID"]:
        data_reco["inner_{0}".format(variable)] = data_reco["{0}".format(variable)]
    
    return data_reco


def get_reco_variables(data_reco, template_type, correction_type):
    if template_type == "Glo" or template_type == "Glo_Trk": # (modify for pt once there are outer corrections)
        for variable in ["pt", "eta", "phi", "charge"]:
            data_reco["reco_{0}".format(variable)] = data_reco["outer_{0}".format(variable)]

    else:
        for variable in ["eta", "phi", "charge"]:
            data_reco["reco_{0}".format(variable)] = data_reco["inner_{0}".format(variable)]
        data_reco["reco_pt"] = data_reco["inner_pt_RC_{0}".format(correction_type)]
    
    return data_reco


def get_paired_outer_value(data_reco): # only for template_type "Sta_Trk"
        # Note: data_reco splitting
        data_reco1 = data_reco[ak.num(data_reco["outer_pt"])==0]
        data_reco2 = data_reco[ak.num(data_reco["outer_pt"])!=0] 
                
         # pairing inner and outer tracks
        for value in ["pt", "eta", "phi", "charge"]:
            pairing_values = ak.cartesian([data_reco2["reco_{0}".format(value)], data_reco2["outer_{0}".format(value)]])
            inner_value, outer_value = ak.unzip(pairing_values)

            counts = ak.flatten(ak.num(outer_value, axis=1) / ak.num(data_reco2["reco_{0}".format(value)], axis=1) * ak.ones_like(data_reco2["reco_{0}".format(value)]))
            counts = ak.values_astype(counts, "int")
            array = outer_value
            data_reco2["paired_outer_{0}".format(value)] = ak.unflatten(array, counts, axis=1)
 
        # define delta_R2 between inner and outer
        data_reco2["delR2_io"] = ((data_reco2["paired_outer_eta"] - data_reco2["reco_eta"])**2 + (data_reco2["paired_outer_phi"] - data_reco2["reco_phi"])**2)
        data_reco2["delR2_io_min"] = ak.min(data_reco2["delR2_io"], axis=-1)

        data_reco1["delR2_io_min"] = ak.ones_like(data_reco1["reco_pt"]) # sets to a value larger than delR2_io cut
        
        # select inner and outer pair
        mask_io_pair = data_reco2["delR2_io_min"] == data_reco2["delR2_io"]
        for value in ["pt", "eta", "phi", "charge"]:
            min_outer_value = data_reco2["paired_outer_{0}".format(value)][mask_io_pair]
            data_reco2["outer_{0}_min".format(value)] = ak.flatten(min_outer_value, axis=-1)
            
            data_reco1["outer_{0}_min".format(value)] = ak.ones_like(data_reco1["reco_pt"]) # designed to fail outer pt in probe mask

        data_reco3 = ak.concatenate((data_reco2, data_reco1))
        
        return data_reco3
    

def get_paired_inner_value(data_reco): # only for template_type "Glo" and "Glo_Trk"
    # Note: data_reco selection
    data_reco = data_reco[ak.num(data_reco["outer_pt"])!=0]
           
    # pairing inner and outer tracks
    for value in ["pt", "eta", "phi", "charge"]:
        pairing_values = ak.cartesian([data_reco["reco_{0}".format(value)], data_reco["inner_{0}".format(value)]])
        outer_value, inner_value = ak.unzip(pairing_values)

        counts = ak.flatten(ak.num(inner_value, axis=1) / ak.num(data_reco["reco_{0}".format(value)], axis=1) * ak.ones_like(data_reco["reco_{0}".format(value)]))
        counts = ak.values_astype(counts, "int")
        array = inner_value
        data_reco["paired_inner_{0}".format(value)] = ak.unflatten(array, counts, axis=1)
         
    # define delta_R2 between inner and outer
    data_reco["delR2_io"] = ((data_reco["paired_inner_eta"] - data_reco["outer_eta"])**2 + (data_reco["paired_inner_phi"] - data_reco["outer_phi"])**2)
    data_reco["delR2_io_min"] = ak.min(data_reco["delR2_io"], axis=-1)
    
    # select inner and outer pair
    mask_io_pair = data_reco["delR2_io_min"] == data_reco["delR2_io"]
    for value in ["pt", "eta", "phi", "charge"]:
        min_inner_value = data_reco["paired_inner_{0}".format(value)][mask_io_pair]
        data_reco["inner_{0}_min".format(value)] = ak.flatten(min_inner_value, axis=-1)
    
    return data_reco


def tag_criteria(data_reco, template_type):
    if template_type == "Glo" or template_type == "Glo_Trk":
        criteria_tag_0 = (data_reco["inner_HLT_trigger_pass"] == True) & (data_reco["inner_Muon_ID"] >= 4)
        
        #reshape (pair criteria tag)
        pair = ak.cartesian([data_reco["outer_pt"], criteria_tag_0])
        outer_value, inner_value = ak.unzip(pair)

        counts = ak.flatten(ak.num(inner_value, axis=1) / ak.num(data_reco["outer_pt"], axis=1) * ak.ones_like(data_reco["outer_pt"]))
        counts = ak.values_astype(counts, "int")
        array = inner_value
        paired_value = ak.unflatten(array, counts, axis=1)
        
        
        mask_io_pair = data_reco["delR2_io_min"] == data_reco["delR2_io"]
        min_inner_value = paired_value[mask_io_pair]
        tag_criteria = ak.flatten(min_inner_value, axis=-1)
        
    else:
        tag_criteria = (data_reco["inner_HLT_trigger_pass"] == True) & (data_reco["inner_Muon_ID"] >= 4)
        
    return tag_criteria


def base_criteria(data_reco, template_type, cut_values):
    if template_type == "HLT":
        base_mask = (data_reco["inner_Muon_ID"] >= 4)
        
    elif template_type == "ID":
        base_mask = (data_reco["inner_Muon_ID"] >= 3)
    
    elif template_type == "Glo":  
        base_mask = (data_reco["inner_pt_min"] > cut_values["pt_inner"]) & (abs(data_reco["inner_eta_min"]) < cut_values["eta_inner"])
    
    elif template_type == "Sta_Trk":
        base_mask = (data_reco["nPixelHits"] > 0 ) & (data_reco["nTrackerLayers"] > 5) & (data_reco["trackAlgo"] != 13) & (data_reco["trackAlgo"] != 14)

    elif template_type == "HLT_Trk":
        track_mask = (data_reco["nPixelHits"] > 0 ) & (data_reco["nTrackerLayers"] > 5) & (data_reco["trackAlgo"] != 13) & (data_reco["trackAlgo"] != 14)
        base_mask = (data_reco["inner_Muon_ID"] >= 4) & track_mask

    elif template_type == "ID_Trk":
        track_mask = (data_reco["nPixelHits"] > 0 ) & (data_reco["nTrackerLayers"] > 5) & (data_reco["trackAlgo"] != 13) & (data_reco["trackAlgo"] != 14)
        base_mask = (data_reco["inner_Muon_ID"] >= 3) & track_mask

    elif template_type == "Glo_Trk":
        track_mask = (data_reco["nPixelHits"] > 0 ) & (data_reco["nTrackerLayers"] > 5) & (data_reco["trackAlgo"] != 13) & (data_reco["trackAlgo"] != 14)
        
        #reshape (pair criteria tag)
        pair = ak.cartesian([data_reco["outer_pt"], track_mask])
        outer_value, inner_value = ak.unzip(pair)

        counts = ak.flatten(ak.num(inner_value, axis=1) / ak.num(data_reco["outer_pt"], axis=1) * ak.ones_like(data_reco["outer_pt"]))
        counts = ak.values_astype(counts, "int")
        array = inner_value
        paired_value = ak.unflatten(array, counts, axis=1)
        
        
        mask_io_pair = data_reco["delR2_io_min"] == data_reco["delR2_io"]
        min_inner_value = paired_value[mask_io_pair]
        reshaped_track_mask = ak.flatten(min_inner_value, axis=-1)
        
        base_mask = (data_reco["inner_pt_min"] > cut_values["pt_inner"]) & (abs(data_reco["inner_eta_min"]) < cut_values["eta_inner"]) & reshaped_track_mask
        
    return base_mask


def probe_criteria_pass(data_reco, template_type, cut_values):
    if template_type == "HLT" or template_type == "HLT_Trk":
        pass_mask = (data_reco["inner_HLT_trigger_pass"] == True)
        
    elif template_type == "ID" or template_type == "ID_Trk":
        pass_mask = (data_reco["inner_Muon_ID"] > 3)
        
    
    elif template_type == "Glo" or template_type == "Glo_Trk": 
        pass_ID_0 = data_reco["inner_Muon_ID"] >= 3
        
        # reshape (pair ID to outer)
        pair = ak.cartesian([data_reco["outer_pt"], pass_ID_0])
        outer_value, inner_value = ak.unzip(pair)

        counts = ak.flatten(ak.num(inner_value, axis=1) / ak.num(data_reco["outer_pt"], axis=1) * ak.ones_like(data_reco["outer_pt"]))
        counts = ak.values_astype(counts, "int")
        array = inner_value
        paired_value = ak.unflatten(array, counts, axis=1)
         
        mask_io_pair = data_reco["delR2_io_min"] == data_reco["delR2_io"]
        min_inner_value = paired_value[mask_io_pair]
        pass_ID = ak.flatten(min_inner_value, axis=-1)
        
        pass_mask = (pass_ID == True) &  (data_reco["delR2_io_min"] < cut_values["delR2_io"])
        
        
    elif template_type == "Sta_Trk":
        pass_mask = (data_reco["delR2_io_min"] < cut_values["delR2_io"]) & (data_reco["outer_pt_min"] > cut_values["pt_outer"]) & (abs(data_reco["outer_eta_min"]) < cut_values["eta_outer"])
    
    return pass_mask


def pair_mask_reco(muon_plus, muon_minus, criteria_1, criteria_2):
    pairs_1 = ak.cartesian([muon_plus[criteria_1], muon_minus[criteria_2]])
    l_idx_1, r_idx_1 = ak.unzip(pairs_1)
    mask_1 = (l_idx_1 == r_idx_1) & (l_idx_1 + r_idx_1 == True)
    
    if criteria_1 == criteria_2: # to prevent double counting
        mask_2 = ak.ones_like(mask_1) == 0
        
    else:
        pairs_2 = ak.cartesian([muon_plus[criteria_2], muon_minus[criteria_1]])
        l_idx_2, r_idx_2 = ak.unzip(pairs_2)
        mask_2 = (l_idx_2 == r_idx_2) & (l_idx_2 + r_idx_2 == True)
    
    return mask_1, mask_2


def pair_mask_gen(data_reco, criteria_1, criteria_2):
    
    l_idx_1, r_idx_1 = data_reco["{0}_plus".format(criteria_1)] & data_reco["gen_cut_plus"], data_reco["{0}_minus".format(criteria_2)] & data_reco["gen_cut_minus"]
    mask_1 = (l_idx_1 == r_idx_1) & (l_idx_1 + r_idx_1 == True)

    if criteria_1 == criteria_2: # to prevent double counting
        mask_2 = ak.ones_like(mask_1) == 0
    
    else:
        l_idx_2, r_idx_2 = data_reco["{0}_plus".format(criteria_2)] & data_reco["gen_cut_plus"], data_reco["{0}_minus".format(criteria_1)] & data_reco["gen_cut_minus"]
        mask_2 = (l_idx_2 == r_idx_2) & (l_idx_2 + r_idx_2 == True)
    
    return mask_1, mask_2


def subtype(template_subtype):
    if template_subtype == "pass":
        criteria_1, criteria_2 = "tag_pass", "probe_pass"
    elif template_subtype == "fail":
        criteria_1, criteria_2 = "tag_pass", "probe_fail"
    elif template_subtype == "2":
        criteria_1, criteria_2 = "probe_pass", "probe_pass"
    elif template_subtype == "1":
        criteria_1, criteria_2 = "probe_pass", "probe_fail"
    elif template_subtype == "0":
        criteria_1, criteria_2 = "probe_fail", "probe_fail"
    
    return criteria_1, criteria_2


def pair_reco_to_gen(data_reco, cut_values, template_type):
    if template_type == "Glo" or template_type == "Glo_Trk":
        cut_eta, cut_pt = cut_values["eta_outer"], cut_values["pt_outer"]
    else:
        cut_eta, cut_pt = cut_values["eta_inner"], cut_values["pt_inner"]
    

    data_reco["reco_cuts"] = (data_reco["reco_eta"] < cut_eta) & (data_reco["reco_pt"] > cut_pt)
    
    data_reco["delR2_minus_min"] = ak.min(data_reco["delR2_minus"], axis=-1)
    data_reco["delR2_plus_min"] = ak.min(data_reco["delR2_plus"], axis=-1)

    # match_gen_reco_minus = data_reco["delR2_minus"] == data_reco["delR2_minus_min"]
    # match_gen_reco_plus = data_reco["delR2_plus"] == data_reco["delR2_plus_min"]
    match_gen_reco_minus = ak.argmin(data_reco["delR2_minus"], axis=-1, keepdims = True)
    match_gen_reco_plus = ak.argmin(data_reco["delR2_plus"], axis=-1, keepdims = True)

    # note change in none value: [False] instead of False if using old match_gen_reco, and axis=0 instead of -1
    for value in ["tag_pass", "probe_pass", "probe_fail", "reco_cuts"]:
        min_minus = data_reco["{0}".format(value)][match_gen_reco_minus]
        min_minus = ak.fill_none(min_minus, value = False, axis=-1)
        data_reco["{0}_minus".format(value)] = ak.flatten(min_minus, axis=-1)

        min_plus = data_reco["{0}".format(value)][match_gen_reco_plus]
        min_plus = ak.fill_none(min_plus, value = False, axis=-1)
        data_reco["{0}_plus".format(value)] = ak.flatten(min_plus, axis=-1)
    
    return data_reco


# arguments for generate_templates:
# template_type = ["HLT", "ID", "Glo", "HLT_Trk", "ID_Trk", "Glo_Trk", "Sta_Trk", "all"]
# template_subtype = ["2", "1", "0", "pass", "fail"]
# correction_type = ["none", "nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]
# cut_type = ["reco", "gen"]
# template_type ending in "_Trk" refers to taking the Muon and Track values for inner values. Obviously, there is only "Sta_Trk" and no "Sta".
# For all "_Trk" types, there are the additional seed and track requirements in the base mask (criteria that need to be satisfied for both pass and fail probes).
# template_type "all" only works with cut_type = "gen", and produces all muons satisfying gen pt and eta cuts. This does not depend on "template_subtype or correction_type, so any value can be entered there.
# template_subtype "2", "1" and "0" refer to 2 passing probes, 1 passing and 1 failing probe, and 2 failing probes respectively. "pass" and "fail" refer to a tag-and-probe with passing and failing probe respectively.
# cut_values.keys() = ["pt_cut", "eta_cut", "delR2_inner", "delR2_outer", "delR2_io", "pt_inner", "pt_outer", "eta_inner", "eta_outer", "upper_pt"]

def generate_templates(data_reco, template_type, template_subtype, correction_type, cut_type, cut_values, bin_in_PV = True):

    # HLT trigger
    data_reco["HLT_trigger_pass"] = data_reco["Muon_triggerBits"]&1 == 1
    
    # prepares inner variables (if template_type ends with "_Trk", concatenate muons and tracks; otherwise just obtains muons)
    if template_type == "HLT_Trk" or template_type == "ID_Trk" or template_type == "Glo_Trk" or template_type == "Sta_Trk":
        data_reco = get_inner_variables_with_trk(data_reco, correction_type)
    else:
        data_reco = get_inner_variables_without_trk(data_reco, correction_type)
    
    
    # obtain reco variables (outer for "Glo" and "Glo_Trk", and inner for all others)
    data_reco = get_reco_variables(data_reco, template_type, correction_type)

    
    # pairing inner and outer variables (only needed for "Glo", "Glo_Trk" and "Sta_Trk")
    # this step obtains delR2_io_min (minimum delR2 between inner and outer), and also the corresponding value of pt, eta, phi and charge
    if template_type == "Sta_Trk": # note there is reordering of data from splitting and recombining data_reco here
        data_reco = get_paired_outer_value(data_reco)         
    elif template_type == "Glo" or template_type == "Glo_Trk": # note the elimination of data with no outer tracks
        data_reco = get_paired_inner_value(data_reco)

    
    if template_type != "all":
        # define base criteria (criteria that has to be true for both pass and fail probes, but not tag)
        base_mask = base_criteria(data_reco, template_type, cut_values)

        # define passing and failing criteria
        probe_pass_mask = probe_criteria_pass(data_reco, template_type, cut_values)
        # probe_fail_mask = probe_criteria_fail(data_reco, template_type)

        data_reco["probe_pass"] = (base_mask == True) & (probe_pass_mask == True)
        data_reco["probe_fail"] = (base_mask == True) & (probe_pass_mask == False)
        # data_reco["probe_fail"] = (base_mask == True) & (probe_fail_mask == True)

        # define tag criteria
        data_reco["tag_pass"] = tag_criteria(data_reco, template_type)


        # reshape nPV into reco-shape
        data_reco["primary_vertices"]  = ak.ones_like(data_reco["reco_pt"]) * data_reco["nPV"]

    
        # delR2 between reco and gen
        data_reco["delR2_minus"] = (data_reco["reco_eta"] - data_reco["lepton_genEta"])**2 + (data_reco["reco_phi"] - data_reco["lepton_genPhi"])**2
        data_reco["delR2_plus"] = (data_reco["reco_eta"] - data_reco["antiLepton_genEta"])**2 + (data_reco["reco_phi"] - data_reco["antiLepton_genPhi"])**2

        # define delR2 cut between reco and gen values
        if template_type == "Glo" or template_type == "Glo_Trk":
            delR2_recogen = cut_values["delR2_outer"]
        else:
            delR2_recogen = cut_values["delR2_inner"]
    
    
    # applyng subtype mask to muon pairs:
    if cut_type == "gen":
        
        if template_type == "all":
            l_idx, r_idx = (abs(data_reco["antiLepton_genEta"]) < cut_values["eta_cut"]) & (data_reco["antiLepton_genPt"] > cut_values["pt_cut"]) & (data_reco["antiLepton_genPt"] < cut_values["upper_pt"]), (abs(data_reco["lepton_genEta"]) < cut_values["eta_cut"]) & (data_reco["lepton_genPt"] > cut_values["pt_cut"]) & (data_reco["lepton_genPt"] < cut_values["upper_pt"])
            mask_1 = (l_idx == r_idx) & (l_idx + r_idx == True)
            mask_2 = ak.ones_like(mask_1) == 0
        
        else:
            # pairing reco to gen variables
            # this step defines delR2_min, and attaches the corresponding mask value of reco
            data_reco = pair_reco_to_gen(data_reco, cut_values, template_type)

            delR2_plus_min_nonone = ak.fill_none(data_reco["delR2_plus_min"], value = 1, axis=0)
            delR2_minus_min_nonone = ak.fill_none(data_reco["delR2_minus_min"], value = 1, axis=0)

            # pt-eta cuts, with delR2 matching requirements between reco and gen
            data_reco["gen_cut_plus"] = (abs(data_reco["antiLepton_genEta"]) < cut_values["eta_cut"]) & (data_reco["antiLepton_genPt"] > cut_values["pt_cut"]) & (data_reco["antiLepton_genPt"] < cut_values["upper_pt"]) & (delR2_plus_min_nonone < delR2_recogen) & (data_reco["reco_cuts_plus"] == True)
            data_reco["gen_cut_minus"] = (abs(data_reco["lepton_genEta"]) < cut_values["eta_cut"]) & (data_reco["lepton_genPt"] > cut_values["pt_cut"]) & (data_reco["lepton_genPt"] < cut_values["upper_pt"]) & (delR2_minus_min_nonone < delR2_recogen) & (data_reco["reco_cuts_minus"] == True)

            # mask for template subtype
            criteria_1, criteria_2 = subtype(template_subtype)
            mask_1, mask_2 = pair_mask_gen(data_reco, criteria_1, criteria_2)
        
        
        # apply subtype mask to plus-minus muon pairs
        mass = np.array([])
        PV = np.array([])
        for mask in (mask_1, mask_2):
            l_eta = data_reco["antiLepton_genEta"][mask]
            r_eta = data_reco["lepton_genEta"][mask]
            l_phi = data_reco["antiLepton_genPhi"][mask]
            r_phi = data_reco["lepton_genPhi"][mask]
            l_pt = data_reco["antiLepton_genPt"][mask]
            r_pt = data_reco["lepton_genPt"][mask]
            l_PV = data_reco["nPV"][mask]
            r_PV = data_reco["nPV"][mask]
            
            muonMass = 0.105658369
            mu1 = vector.arr({"pt":l_pt, "phi":l_phi, "eta": l_eta, "mass": l_pt*0+muonMass})
            mu2 = vector.arr({"pt":r_pt, "phi":r_phi, "eta": r_eta, "mass": r_pt*0+muonMass})

            mass = np.concatenate((mass, (mu1+mu2).mass))
            PV = np.concatenate((PV, l_PV))
        
    
    elif cut_type == "reco":
        # split into plus and minus muons, along with pt-eta cuts and delR2 matching requirements
        mask_plus  = (data_reco["delR2_plus"] < delR2_recogen) & (data_reco["reco_charge"] == 1)
        mask_minus  = (data_reco["delR2_minus"] < delR2_recogen) & (data_reco["reco_charge"] == -1)
        cut_mask = (abs(data_reco["reco_eta"]) < cut_values["eta_cut"]) & (data_reco["reco_pt"] > cut_values["pt_cut"]) & (data_reco["reco_pt"] < cut_values["upper_pt"])

        muon_plus = data_reco[["tag_pass", "probe_pass", "probe_fail", "reco_pt", "reco_eta", "reco_phi", "primary_vertices"]][mask_plus & cut_mask]
        muon_minus = data_reco[["tag_pass", "probe_pass", "probe_fail", "reco_pt", "reco_eta", "reco_phi", "primary_vertices"]][mask_minus & cut_mask]

    
        # mask for template subtype
        criteria_1, criteria_2 = subtype(template_subtype)
        mask_1, mask_2 = pair_mask_reco(muon_plus, muon_minus, criteria_1, criteria_2)   
        
        # apply subtype mask to plus-minus muon pairs
        mass = np.array([])
        PV = np.array([])
        for mask in (mask_1, mask_2):

            pairs_eta = ak.cartesian([muon_plus["reco_eta"], muon_minus["reco_eta"]])[mask]
            pairs_phi = ak.cartesian([muon_plus["reco_phi"], muon_minus["reco_phi"]])[mask]
            pairs_pt = ak.cartesian([muon_plus["reco_pt"], muon_minus["reco_pt"]])[mask]
            pairs_PV = ak.cartesian([muon_plus["primary_vertices"], muon_minus["primary_vertices"]])[mask]

            l_idx_eta, r_idx_eta = ak.unzip(pairs_eta)
            l_idx_phi, r_idx_phi = ak.unzip(pairs_phi)
            l_idx_pt, r_idx_pt = ak.unzip(pairs_pt)
            l_idx_PV, r_idx_PV = ak.unzip(pairs_PV)

            l_eta = ak.flatten(l_idx_eta)
            r_eta = ak.flatten(r_idx_eta)
            l_phi = ak.flatten(l_idx_phi)
            r_phi = ak.flatten(r_idx_phi)
            l_pt = ak.flatten(l_idx_pt)
            r_pt = ak.flatten(r_idx_pt)
            l_PV = ak.flatten(l_idx_PV)
            r_PV = ak.flatten(r_idx_PV)

            muonMass = 0.105658369
            mu1 = vector.arr({"pt":l_pt, "phi":l_phi, "eta": l_eta, "mass": l_pt*0+muonMass})
            mu2 = vector.arr({"pt":r_pt, "phi":r_phi, "eta": r_eta, "mass": r_pt*0+muonMass})


            mass = np.concatenate((mass, (mu1+mu2).mass))
            PV = np.concatenate((PV, l_PV))

            
    # binning into mass (and PV if bin_in_PV == True (default))
    if template_type == "Glo" or template_type == "Glo_Trk":
        nbins, bin_min, bin_max = 70, 56, 126
    else:
        nbins, bin_min, bin_max = 50, 66, 116
        
    if bin_in_PV == True:
        hist = Hist.new.Regular(nbins, bin_min, bin_max, name="mass").Regular(100, 0.5, 100.5, name="PV").Double()
        hist.fill(mass, PV)
    else:
        hist = Hist.new.Regular(nbins, bin_min, bin_max, name="mass").Double()
        hist.fill(mass)    
     
    return hist.view()


#######################

parser = argparse.ArgumentParser()
parser.add_argument("mc_data_file")
parser.add_argument("rochester_file")
args = parser.parse_args()


filename_ntuples = args.mc_data_file
rochester_filename = args.rochester_file


f1 = uproot.open(filename_ntuples)

treename = "zcounting/tree"
data = uproot.open(filename_ntuples+":"+treename)


columns = ["lepton_genPt","lepton_genEta", "lepton_genPhi",
 "antiLepton_genPt","antiLepton_genEta","antiLepton_genPhi",
 "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge",
 "Muon_ID", "Muon_triggerBits", "nPV",
 "Muon_useUpdated", "Muon_ptStaUpd", "Muon_etaStaUpd", "Muon_phiStaUpd", "Muon_chargeStaUpd",
 "Muon_ptStaReg", "Muon_etaStaReg", "Muon_phiStaReg", "Muon_chargeStaReg",
 "Track_pt", "Track_eta", "Track_phi", "Track_charge",
 "Muon_ptTrk", "Muon_etaTrk", "Muon_phiTrk", "Muon_chargeTrk",
 "Track_nPixelHits", "Track_nTrackerLayers", "Track_trackAlgo",
 'Muon_nPixelHits', 'Muon_nTrackerLayers','Muon_trackAlgo', "decayMode", "nMuon"]


keys_reco1 = []
for correction in ["nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]:
    _add_keys_reco1 = ["HLT2_{0}".format(correction), "HLT1_{0}".format(correction), "IDfail_{0}".format(correction), "Stapass_{0}".format(correction), "Stafail_{0}".format(correction)]
    keys_reco1 = keys_reco1 + _add_keys_reco1
keys_reco2 = ["Glopass_none", "Glofail_none"]
keys_gen1 = ["HLT2", "HLT1", "HLT0", "ID2", "ID1", "ID0", "IDpass", "IDfail",
 "Sta_Trk2", "Sta_Trk1", "Sta_Trk0", "Sta_Trkpass", "Sta_Trkfail", "HLT_Trkpass", "HLT_Trkfail", "ID_Trkpass", "ID_Trkfail", "all"]
keys_gen2 = ["Glo2", "Glo1", "Glo0", "Glopass", "Glofail", "Glo_Trkfail"]

rochester = _get_rochester_corr(rochester_filename)


result_gen = {}
result_reco = {}
for key in keys_reco1:
    result_reco[key] = np.zeros((50,100))
for key in keys_reco2:
    result_reco[key] = np.zeros((70,100))
for key in keys_gen1:
    result_gen[key] = np.zeros((50,100))
for key in keys_gen2:
    result_gen[key] = np.zeros((70,100))

keys_gen = keys_gen1 + keys_gen2
keys_reco = keys_reco1 + keys_reco2

cut_values = {"pt_cut": 27, "eta_cut": 2.4, "delR2_inner": 0.0025, "delR2_outer": 0.09, "delR2_io": 0.09, "pt_inner": 15, "pt_outer": 15, "eta_inner": 2.5, "eta_outer": 2.4, "upper_pt":200}


def produce_all_templates(data_reco, rochester, cut_values, keys_gen, keys_reco, result_gen, result_reco):
    
    # define outer pt, eta, phi, charge
    data_reco["outer_pt"] = data_reco["Muon_useUpdated"] * data_reco["Muon_ptStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_ptStaReg"] 
    data_reco["outer_eta"] = data_reco["Muon_useUpdated"] * data_reco["Muon_etaStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_etaStaReg"] 
    data_reco["outer_phi"] = data_reco["Muon_useUpdated"] * data_reco["Muon_phiStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_phiStaReg"] 
    data_reco["outer_charge"] = data_reco["Muon_useUpdated"] * data_reco["Muon_chargeStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_chargeStaReg"]

    # get corrections
    data_reco1 = get_corrections(data_reco, rochester)

    # produce templates
    #reco templates
    reco_templates = {}
    for template_type in ["HLT"]:
        for template_subtype in ["2", "1"]:
            for correction_type in ["nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]:
                reco_templates["{0}{1}_{2}".format(template_type, template_subtype, correction_type)] = generate_templates(data_reco1, template_type, template_subtype, correction_type, "reco", cut_values)

    for template_type in ["ID"]:
        for template_subtype in ["fail"]:
            for correction_type in ["nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]:
                reco_templates["{0}{1}_{2}".format(template_type, template_subtype, correction_type)] = generate_templates(data_reco1, template_type, template_subtype, correction_type, "reco", cut_values)
        
    for template_type in ["Glo"]: # to include correction_types once they are obtained
        for template_subtype in ["pass", "fail"]:
            reco_templates["{0}{1}_none".format(template_type, template_subtype)] = generate_templates(data_reco1, template_type, template_subtype, "none", "reco", cut_values)
    
    for template_type in ["Sta_Trk"]:
        for template_subtype in ["pass", "fail"]:
            for correction_type in ["nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]:
                reco_templates["Sta{1}_{2}".format(template_type, template_subtype, correction_type)] = generate_templates(data_reco1, template_type, template_subtype, correction_type, "reco", cut_values)
    
    #gen templates
    gen_templates = {}
    for template_type in ["HLT"]:
        for template_subtype in ["2", "1", "0"]:
            gen_templates["{0}{1}".format(template_type, template_subtype)] = generate_templates(data_reco1, template_type, template_subtype, "nominal", "gen", cut_values)

    for template_type in ["ID", "Sta_Trk"]:
        for template_subtype in ["pass", "fail", "2", "1", "0"]:
            gen_templates["{0}{1}".format(template_type, template_subtype)] = generate_templates(data_reco1, template_type, template_subtype, "nominal", "gen", cut_values)
        
    for template_type in ["Glo"]:
        for template_subtype in ["pass", "fail", "2", "1", "0"]:
            gen_templates["{0}{1}".format(template_type, template_subtype)] = generate_templates(data_reco1, template_type, template_subtype, "none", "gen", cut_values)

    for template_type in ["HLT_Trk", "ID_Trk"]:
        for template_subtype in ["pass", "fail"]:
            gen_templates["{0}{1}".format(template_type, template_subtype)] = generate_templates(data_reco1, template_type, template_subtype, "nominal", "gen", cut_values)

    for template_type in ["Glo_Trk"]:
        for template_subtype in ["fail"]:
            gen_templates["{0}{1}".format(template_type, template_subtype)] = generate_templates(data_reco1, template_type, template_subtype, "none", "gen", cut_values)
    
    gen_templates["all"] = generate_templates(data_reco1, "all", "pass", "none", "gen", cut_values)
    
    
    # add into total templates
    for key in keys_gen:
        result_gen[key] = result_gen[key] + gen_templates[key]
    
    for key in keys_reco:
        result_reco[key] = result_reco[key] + reco_templates[key]

    return result_gen, result_reco

i = 0
for array in uproot.iterate(filename_ntuples+":"+treename, expressions = columns, step_size=100000):
    data_reco = array[columns]
    data_reco = data_reco[(data_reco["decayMode"]==13)]

    result_gen, result_reco = produce_all_templates(data_reco, rochester, cut_values, keys_gen, keys_reco, result_gen, result_reco)
    
    i = i+1
    print(i)


with h5py.File('templates/templates_gen.hdf5', 'w') as outfile:
    for dset_name in result_gen:
        dset = outfile.create_dataset(dset_name, data = result_gen[dset_name])

with h5py.File('templates/templates_reco.hdf5', 'w') as outfile:
    for dset_name in result_reco:
        dset = outfile.create_dataset(dset_name, data = result_reco[dset_name])

