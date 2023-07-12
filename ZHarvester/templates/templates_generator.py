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


#####################################

def generate_templates_HLT(data_reco, rochester, bin_in_PV = True):
      
        
    # define delta_R2
    data_reco["delR2_minus"] = (data_reco["Muon_eta"] - data_reco["muon_genEta"])**2 + (data_reco["Muon_phi"] - data_reco["muon_genPhi"])**2
    data_reco["delR2_plus"] = (data_reco["Muon_eta"] - data_reco["antiMuon_genEta"])**2 + (data_reco["Muon_phi"] - data_reco["antiMuon_genPhi"])**2

    # format nPV into usable shape for masks
    data_reco["primary_vertices"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["nPV"]
    
    
    data_reco["correct_charge_genPt"] = (data_reco["Muon_charge"] + 1)/2 * data_reco["antiMuon_genPt"] - (data_reco["Muon_charge"] - 1)/2 * data_reco["muon_genPt"]

    # format gen variables into usable shape for masks
    data_reco["gen_pt"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["muon_genPt"]
    data_reco["gen_eta"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["muon_genEta"]
    data_reco["gen_phi"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["muon_genPhi"]
    data_reco["anti_gen_pt"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["antiMuon_genPt"]
    data_reco["anti_gen_eta"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["antiMuon_genEta"]
    data_reco["anti_gen_phi"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["antiMuon_genPhi"]

    
    # define masks, except pt requirement
    mask_minus = (data_reco["delR2_minus"]<0.01) & (data_reco["Muon_charge"] == -1) & (abs(data_reco["Muon_eta"]) < 2.4) & (data_reco["Muon_ID"] >= 4) 
    mask_plus  = (data_reco["delR2_plus"]<0.01)  & (data_reco["Muon_charge"] == 1)  & (abs(data_reco["Muon_eta"]) < 2.4) & (data_reco["Muon_ID"] >= 4) 
    
    
    data_reco["HLT_trigger_pass"] = data_reco["Muon_triggerBits"]&1 == 1
    
    
    # apply mask for all reconstructed muons:
    muon_plus = data_reco[["Muon_charge", "Muon_eta","Muon_phi", "Muon_pt", "HLT_trigger_pass", "primary_vertices", "correct_charge_genPt", "anti_gen_eta", "anti_gen_phi", "anti_gen_pt"]][mask_plus]
    muon_minus = data_reco[["Muon_charge", "Muon_eta","Muon_phi", "Muon_pt", "HLT_trigger_pass", "primary_vertices", "correct_charge_genPt", "gen_eta", "gen_phi", "gen_pt"]][mask_minus]
        
     
        
    # apply rochester correction, and its errors
        
    scale_plus = rochester.kScaleDT(muon_plus["Muon_charge"], muon_plus["Muon_pt"], muon_plus["Muon_eta"], muon_plus["Muon_phi"])
    scale_minus = rochester.kScaleDT(muon_minus["Muon_charge"], muon_minus["Muon_pt"], muon_minus["Muon_eta"], muon_minus["Muon_phi"])
    spread_plus = rochester.kSpreadMC(muon_plus["Muon_charge"], muon_plus["Muon_pt"], muon_plus["Muon_eta"], muon_plus["Muon_phi"], muon_plus["correct_charge_genPt"])
    spread_minus = rochester.kSpreadMC(muon_minus["Muon_charge"], muon_minus["Muon_pt"], muon_minus["Muon_eta"], muon_minus["Muon_phi"], muon_minus["correct_charge_genPt"])

    scale_plus_error = rochester.kScaleDTerror(muon_plus["Muon_charge"], muon_plus["Muon_pt"], muon_plus["Muon_eta"], muon_plus["Muon_phi"])
    scale_minus_error = rochester.kScaleDTerror(muon_minus["Muon_charge"], muon_minus["Muon_pt"], muon_minus["Muon_eta"], muon_minus["Muon_phi"])
    spread_plus_error = rochester.kSpreadMCerror(muon_plus["Muon_charge"], muon_plus["Muon_pt"], muon_plus["Muon_eta"], muon_plus["Muon_phi"], muon_plus["correct_charge_genPt"])
    spread_minus_error = rochester.kSpreadMCerror(muon_minus["Muon_charge"], muon_minus["Muon_pt"], muon_minus["Muon_eta"], muon_minus["Muon_phi"], muon_minus["correct_charge_genPt"])

   
    # all the templates
    template_types = ["none", "nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]
    result = {}


    muon_plus["Muon_pt_RC_none"] = muon_plus["Muon_pt"]
    muon_minus["Muon_pt_RC_none"] = muon_minus["Muon_pt"]
    
    muon_plus["Muon_pt_RC_nominal"] = muon_plus["Muon_pt"] * scale_plus * spread_plus
    muon_minus["Muon_pt_RC_nominal"] = muon_minus["Muon_pt"] * scale_minus * spread_minus

    muon_plus["Muon_pt_RC_plus_scale"] = muon_plus["Muon_pt"] * (scale_plus + scale_plus_error) * spread_plus
    muon_minus["Muon_pt_RC_plus_scale"] = muon_minus["Muon_pt"] * (scale_minus + scale_minus_error) * spread_minus

    muon_plus["Muon_pt_RC_minus_scale"] = muon_plus["Muon_pt"] * (scale_plus - scale_plus_error) * spread_plus
    muon_minus["Muon_pt_RC_minus_scale"] = muon_minus["Muon_pt"] * (scale_minus - scale_minus_error) * spread_minus

    muon_plus["Muon_pt_RC_plus_spread"] = muon_plus["Muon_pt"] * scale_plus * (spread_plus + spread_plus_error)
    muon_minus["Muon_pt_RC_plus_spread"] = muon_minus["Muon_pt"] * scale_minus * (spread_minus + spread_minus_error)

    muon_plus["Muon_pt_RC_minus_spread"] = muon_plus["Muon_pt"] * scale_plus * (spread_plus - spread_plus_error)
    muon_minus["Muon_pt_RC_minus_spread"] = muon_minus["Muon_pt"] * scale_minus * (spread_minus - spread_minus_error)


    for rochester_error_type in template_types:
            
        # mask for pt
        muon_reco_plus = muon_plus[muon_plus["Muon_pt_RC_{0}".format(rochester_error_type)] > 25]
        muon_reco_minus = muon_minus[muon_minus["Muon_pt_RC_{0}".format(rochester_error_type)] > 25]
            
        
        
        # define masks for HLT:

        pairs_HLT = ak.cartesian([muon_reco_plus["HLT_trigger_pass"], muon_reco_minus["HLT_trigger_pass"]])
        l_idx_HLT, r_idx_HLT = ak.unzip(pairs_HLT)

        mask_HLT2 = (l_idx_HLT == r_idx_HLT) & (l_idx_HLT + r_idx_HLT == True)
        mask_HLT0 = (l_idx_HLT + r_idx_HLT)==False
        mask_HLT1 = (l_idx_HLT != r_idx_HLT)
        no_mask_HLT = ak.ones_like(l_idx_HLT)==1

        HLT_trigger_types = ["HLT2", "HLT1", "HLT0", "noHLT"]

        for HLT_trigger in HLT_trigger_types:

            if HLT_trigger == "HLT2":
                mask_HLT = mask_HLT2
            elif HLT_trigger == "HLT1":
                mask_HLT = mask_HLT1
            elif HLT_trigger == "HLT0":
                mask_HLT = mask_HLT0
            elif HLT_trigger == "noHLT":
                mask_HLT = no_mask_HLT
                
            cuts = ["gen_cut", "reco_cut"]
            for cut in cuts:

                if cut == "gen_cut":
                    eta_values_plus, eta_values_minus, phi_values_plus, phi_values_minus, pt_values_plus, pt_values_minus = "anti_gen_eta", "gen_eta", "anti_gen_phi", "gen_phi", "anti_gen_pt", "gen_pt"
                elif cut == "reco_cut":
                    eta_values_plus, eta_values_minus, phi_values_plus, phi_values_minus, pt_values_plus, pt_values_minus = "Muon_eta", 'Muon_eta', "Muon_phi", "Muon_phi", "Muon_pt_RC_{0}".format(rochester_error_type), "Muon_pt_RC_{0}".format(rochester_error_type)


                # pairing, and applying HLT masks:
                pairs_eta = ak.cartesian([muon_reco_plus[eta_values_plus], muon_reco_minus[eta_values_minus]])[mask_HLT]
                pairs_phi = ak.cartesian([muon_reco_plus[phi_values_plus], muon_reco_minus[phi_values_minus]])[mask_HLT]
                pairs_pt = ak.cartesian([muon_reco_plus[pt_values_plus], muon_reco_minus[pt_values_minus]])[mask_HLT]
                pairs_PV = ak.cartesian([muon_reco_plus["primary_vertices"], muon_reco_minus["primary_vertices"]])[mask_HLT]

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

                mass = (mu1+mu2).mass

                PV = l_PV


                # binning
                if bin_in_PV == True:
                    hist = Hist.new.Regular(50, 66, 116, name="mass").Regular(100, 0.5, 100.5, name="PV").Double()
                    hist.fill(mass, PV)
                    result["hist_{0}_{1}_{2}".format(rochester_error_type, HLT_trigger, cut)] = hist.view()
                else:
                    hist = Hist.new.Regular(50, 66, 116, name="mass").Double()
                    hist.fill(mass)
                    result["hist_{0}_{1}_{2}".format(rochester_error_type, HLT_trigger, cut)] =  hist.view()           

    
    rochester_corrections = {}
    rochester_corrections["scale_plus"] = scale_plus
    rochester_corrections["scale_minus"] = scale_minus
    rochester_corrections["spread_plus"] = spread_plus
    rochester_corrections["spread_minus"] = spread_minus
    rochester_corrections["scale_plus_error"] = scale_plus_error
    rochester_corrections["scale_minus_error"] = scale_minus_error
    rochester_corrections["spread_plus_error"] = spread_plus_error
    rochester_corrections["spread_minus_error"] = spread_minus_error
    
    return result, rochester_corrections


#####################################

def generate_templates_ID(data_reco, rochester, bin_in_PV = True):
      
        
    # define delta_R2
    data_reco["delR2_minus"] = (data_reco["Muon_eta"] - data_reco["muon_genEta"])**2 + (data_reco["Muon_phi"] - data_reco["muon_genPhi"])**2
    data_reco["delR2_plus"] = (data_reco["Muon_eta"] - data_reco["antiMuon_genEta"])**2 + (data_reco["Muon_phi"] - data_reco["antiMuon_genPhi"])**2

    # format nPV into usable shape for masks
    data_reco["primary_vertices"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["nPV"]
    
    
    data_reco["correct_charge_genPt"] = (data_reco["Muon_charge"] + 1)/2 * data_reco["antiMuon_genPt"] - (data_reco["Muon_charge"] - 1)/2 * data_reco["muon_genPt"]

    
    # format gen variables into usable shape for masks
    data_reco["gen_pt"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["muon_genPt"]
    data_reco["gen_eta"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["muon_genEta"]
    data_reco["gen_phi"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["muon_genPhi"]
    data_reco["anti_gen_pt"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["antiMuon_genPt"]
    data_reco["anti_gen_eta"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["antiMuon_genEta"]
    data_reco["anti_gen_phi"]  = ak.ones_like(data_reco["Muon_eta"]) * data_reco["antiMuon_genPhi"]

    
    # define masks, except pt requirement
    mask_minus = (data_reco["delR2_minus"]<0.01) & (data_reco["Muon_charge"] == -1) & (abs(data_reco["Muon_eta"]) < 2.4)
    mask_plus  = (data_reco["delR2_plus"]<0.01)  & (data_reco["Muon_charge"] == 1)  & (abs(data_reco["Muon_eta"]) < 2.4)

     
    data_reco["HLT_trigger_pass"] = data_reco["Muon_triggerBits"]&1 == 1
    data_reco["criteria_tag"] = (data_reco["HLT_trigger_pass"] == True) & (data_reco["Muon_ID"] >= 4)
    data_reco["criteria_probe_fail"] = data_reco["Muon_ID"] == 3

    data_reco["ID_pass"] = data_reco["Muon_ID"] >= 4
    
    
    # apply mask for all reconstructed muons:
    muon_plus = data_reco[["Muon_charge", "Muon_eta","Muon_phi", "Muon_pt", "HLT_trigger_pass", "primary_vertices", "correct_charge_genPt", "anti_gen_eta", "anti_gen_phi", "anti_gen_pt", "criteria_tag", "criteria_probe_fail", "ID_pass"]][mask_plus]
    muon_minus = data_reco[["Muon_charge", "Muon_eta","Muon_phi", "Muon_pt", "HLT_trigger_pass", "primary_vertices", "correct_charge_genPt", "gen_eta", "gen_phi", "gen_pt", "criteria_tag", "criteria_probe_fail", "ID_pass"]][mask_minus]    
     
        
    # apply rochester correction, and its errors
        
    scale_plus = rochester.kScaleDT(muon_plus["Muon_charge"], muon_plus["Muon_pt"], muon_plus["Muon_eta"], muon_plus["Muon_phi"])
    scale_minus = rochester.kScaleDT(muon_minus["Muon_charge"], muon_minus["Muon_pt"], muon_minus["Muon_eta"], muon_minus["Muon_phi"])
    spread_plus = rochester.kSpreadMC(muon_plus["Muon_charge"], muon_plus["Muon_pt"], muon_plus["Muon_eta"], muon_plus["Muon_phi"], muon_plus["correct_charge_genPt"])
    spread_minus = rochester.kSpreadMC(muon_minus["Muon_charge"], muon_minus["Muon_pt"], muon_minus["Muon_eta"], muon_minus["Muon_phi"], muon_minus["correct_charge_genPt"])

    scale_plus_error = rochester.kScaleDTerror(muon_plus["Muon_charge"], muon_plus["Muon_pt"], muon_plus["Muon_eta"], muon_plus["Muon_phi"])
    scale_minus_error = rochester.kScaleDTerror(muon_minus["Muon_charge"], muon_minus["Muon_pt"], muon_minus["Muon_eta"], muon_minus["Muon_phi"])
    spread_plus_error = rochester.kSpreadMCerror(muon_plus["Muon_charge"], muon_plus["Muon_pt"], muon_plus["Muon_eta"], muon_plus["Muon_phi"], muon_plus["correct_charge_genPt"])
    spread_minus_error = rochester.kSpreadMCerror(muon_minus["Muon_charge"], muon_minus["Muon_pt"], muon_minus["Muon_eta"], muon_minus["Muon_phi"], muon_minus["correct_charge_genPt"])

   
    # all the templates
    template_types = ["none", "nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]
    result = {}


    muon_plus["Muon_pt_RC_none"] = muon_plus["Muon_pt"]
    muon_minus["Muon_pt_RC_none"] = muon_minus["Muon_pt"]
    
    muon_plus["Muon_pt_RC_nominal"] = muon_plus["Muon_pt"] * scale_plus * spread_plus
    muon_minus["Muon_pt_RC_nominal"] = muon_minus["Muon_pt"] * scale_minus * spread_minus

    muon_plus["Muon_pt_RC_plus_scale"] = muon_plus["Muon_pt"] * (scale_plus + scale_plus_error) * spread_plus
    muon_minus["Muon_pt_RC_plus_scale"] = muon_minus["Muon_pt"] * (scale_minus + scale_minus_error) * spread_minus

    muon_plus["Muon_pt_RC_minus_scale"] = muon_plus["Muon_pt"] * (scale_plus - scale_plus_error) * spread_plus
    muon_minus["Muon_pt_RC_minus_scale"] = muon_minus["Muon_pt"] * (scale_minus - scale_minus_error) * spread_minus

    muon_plus["Muon_pt_RC_plus_spread"] = muon_plus["Muon_pt"] * scale_plus * (spread_plus + spread_plus_error)
    muon_minus["Muon_pt_RC_plus_spread"] = muon_minus["Muon_pt"] * scale_minus * (spread_minus + spread_minus_error)

    muon_plus["Muon_pt_RC_minus_spread"] = muon_plus["Muon_pt"] * scale_plus * (spread_plus - spread_plus_error)
    muon_minus["Muon_pt_RC_minus_spread"] = muon_minus["Muon_pt"] * scale_minus * (spread_minus - spread_minus_error)


    for rochester_error_type in template_types:
            
        # mask for pt
        muon_glo_plus = muon_plus[muon_plus["Muon_pt_RC_{0}".format(rochester_error_type)] > 25]
        muon_glo_minus = muon_minus[muon_minus["Muon_pt_RC_{0}".format(rochester_error_type)] > 25]
            
        
        # define masks for IDfail:

        pairs_1 = ak.cartesian([muon_glo_plus["criteria_tag"], muon_glo_minus["criteria_probe_fail"]])
        l_idx_1, r_idx_1 = ak.unzip(pairs_1)
        mask_1 = (l_idx_1 == r_idx_1) & (l_idx_1 + r_idx_1 == True)

        pairs_2 = ak.cartesian([muon_glo_plus["criteria_probe_fail"], muon_glo_minus["criteria_tag"]])
        l_idx_2, r_idx_2 = ak.unzip(pairs_2)
        mask_2 = (l_idx_2 == r_idx_2) & (l_idx_2 + r_idx_2 == True)

        mask_IDfail = mask_1 + mask_2

        # define masks for ID0, ID1, ID2:
        pairs_ID = ak.cartesian([muon_glo_plus["ID_pass"], muon_glo_minus["ID_pass"]])
        l_idx_ID, r_idx_ID = ak.unzip(pairs_ID)

        mask_ID2 = (l_idx_ID == r_idx_ID) & (l_idx_ID + r_idx_ID == True)
        mask_ID0 = (l_idx_ID + r_idx_ID)==False
        mask_ID1 = (l_idx_ID != r_idx_ID)
        no_mask_ID = ak.ones_like(l_idx_ID)==1


        ID_mask_types = ["ID2", "ID1", "ID0", "IDfail"]

        for ID_type in ID_mask_types:

            if ID_type == "ID2":
                mask_ID = mask_ID2
            elif ID_type == "ID1":
                mask_ID = mask_ID1
            elif ID_type == "ID0":
                mask_ID = mask_ID0
            elif ID_type == "IDfail":
                mask_ID = mask_IDfail       

            cuts = ["gen_cut", "reco_cut"]
            for cut in cuts:

                if cut == "gen_cut":
                    eta_values_plus, eta_values_minus, phi_values_plus, phi_values_minus, pt_values_plus, pt_values_minus = "anti_gen_eta", "gen_eta", "anti_gen_phi", "gen_phi", "anti_gen_pt", "gen_pt"
                elif cut == "reco_cut":
                    eta_values_plus, eta_values_minus, phi_values_plus, phi_values_minus, pt_values_plus, pt_values_minus = "Muon_eta", 'Muon_eta', "Muon_phi", "Muon_phi", "Muon_pt_RC_{0}".format(rochester_error_type), "Muon_pt_RC_{0}".format(rochester_error_type)


                # pairing, and applying HLT masks:
                pairs_eta = ak.cartesian([muon_glo_plus[eta_values_plus], muon_glo_minus[eta_values_minus]])[mask_ID]
                pairs_phi = ak.cartesian([muon_glo_plus[phi_values_plus], muon_glo_minus[phi_values_minus]])[mask_ID]
                pairs_pt = ak.cartesian([muon_glo_plus[pt_values_plus], muon_glo_minus[pt_values_minus]])[mask_ID]
                pairs_PV = ak.cartesian([muon_glo_plus["primary_vertices"], muon_glo_minus["primary_vertices"]])[mask_ID]

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

                mass = (mu1+mu2).mass

                PV = l_PV
        

                # binning
                if bin_in_PV == True:
                    hist = Hist.new.Regular(50, 66, 116, name="mass").Regular(100, 0.5, 100.5, name="PV").Double()
                    hist.fill(mass, PV)
                    result["hist_{0}_{1}_{2}".format(rochester_error_type, ID_type, cut)] = hist.view()
                else:
                    hist = Hist.new.Regular(50, 66, 116, name="mass").Double()
                    hist.fill(mass)
                    result["hist_{0}_{1}_{2}".format(rochester_error_type, ID_type, cut)] =  hist.view()           

    
    rochester_corrections = {}
    rochester_corrections["scale_plus"] = scale_plus
    rochester_corrections["scale_minus"] = scale_minus
    rochester_corrections["spread_plus"] = spread_plus
    rochester_corrections["spread_minus"] = spread_minus
    rochester_corrections["scale_plus_error"] = scale_plus_error
    rochester_corrections["scale_minus_error"] = scale_minus_error
    rochester_corrections["spread_plus_error"] = spread_plus_error
    rochester_corrections["spread_minus_error"] = spread_minus_error
    
    return result, rochester_corrections


#####################################

def generate_templates_Glo(data_reco, rochester, bin_in_PV = True): # note binning resolution

    # define outer pt, eta, phi, charge
    data_reco["outer_pt"] = data_reco["Muon_useUpdated"] * data_reco["Muon_ptStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_ptStaReg"] 
    data_reco["outer_eta"] = data_reco["Muon_useUpdated"] * data_reco["Muon_etaStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_etaStaReg"] 
    data_reco["outer_phi"] = data_reco["Muon_useUpdated"] * data_reco["Muon_phiStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_phiStaReg"] 
    data_reco["outer_charge"] = data_reco["Muon_useUpdated"] * data_reco["Muon_chargeStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_chargeStaReg"] 
      

    data_reco2 = {}    
    # pairing inner and outer tracks
    kinematics = ["pt", "eta", "phi", "charge"]
    for value in kinematics:
        data_reco2["pairing_{0}".format(value)] = ak.cartesian([data_reco["Muon_{0}".format(value)], data_reco["outer_{0}".format(value)]])
        data_reco["paired_inner_{0}".format(value)], data_reco["paired_outer_{0}".format(value)] = ak.unzip(data_reco2["pairing_{0}".format(value)])
    
    # define delta_R2 between inner and outer
    data_reco["delR2_io"] = (data_reco["paired_outer_eta"] - data_reco["paired_inner_eta"])**2 + (data_reco["paired_outer_phi"] - data_reco["paired_inner_phi"])**2

    
    # define tag criteria (reshaped further below)
    data_reco["HLT_trigger_pass"] = data_reco["Muon_triggerBits"]&1 == 1
    data_reco["criteria_tag"] = (data_reco["HLT_trigger_pass"] == True) & (data_reco["Muon_ID"] >= 4)

          
    # reshape gen-type variables
    variables_type_1 = ["nPV", "muon_genPt", "muon_genEta", "muon_genPhi", "antiMuon_genPt", "antiMuon_genEta", "antiMuon_genPhi"]
    for variable in variables_type_1:
        data_reco["paired_{0}".format(variable)]  = ak.ones_like(data_reco["delR2_io"]) * data_reco["{0}".format(variable)]     
    
    # reshape inner-type variables
    variables_type_2 = ["criteria_tag", "Muon_ID"]
    for variable in variables_type_2:
        data_reco2["pairing_{0}".format(variable)] = ak.cartesian([data_reco["{0}".format(variable)], data_reco["outer_pt".format(variable)]])
        data_reco["paired_{0}".format(variable)], data_reco2["paired_{0}".format(variable)] = ak.unzip(data_reco2["pairing_{0}".format(variable)])

    
    # define delta_R2 with gen
    data_reco["delR2_minus"] = (data_reco["paired_outer_eta"] - data_reco["muon_genEta"])**2 + (data_reco["paired_outer_phi"] - data_reco["muon_genPhi"])**2
    data_reco["delR2_plus"] = (data_reco["paired_outer_eta"] - data_reco["antiMuon_genEta"])**2 + (data_reco["paired_outer_phi"] - data_reco["antiMuon_genPhi"])**2
        
    
    data_reco["correct_charge_genPt"] = (data_reco["paired_inner_charge"] + 1)/2 * data_reco["antiMuon_genPt"] - (data_reco["paired_inner_charge"] - 1)/2 * data_reco["muon_genPt"]

    
    # rochester correction for inner pt
    scale = rochester.kScaleDT(data_reco["paired_inner_charge"], data_reco["paired_inner_pt"], data_reco["paired_inner_eta"], data_reco["paired_inner_phi"])
    spread = rochester.kSpreadMC(data_reco["paired_inner_charge"], data_reco["paired_inner_pt"], data_reco["paired_inner_eta"], data_reco["paired_inner_phi"], data_reco["correct_charge_genPt"])

    data_reco["Muon_pt_RC"] = data_reco["paired_inner_pt"] * scale * spread
    
    
    # define masks, except outer pt requirement
    mask_minus = (data_reco["delR2_minus"]<0.09) & (data_reco["delR2_io"]<0.09) & (data_reco["paired_outer_charge"] == -1) & (abs(data_reco["paired_outer_eta"]) < 2.4) & (data_reco["Muon_pt_RC"] > 20) & (abs(data_reco["paired_inner_eta"]) < 2.5) 
    mask_plus  = (data_reco["delR2_plus"]<0.09)  & (data_reco["delR2_io"]<0.09) & (data_reco["paired_outer_charge"] == 1)  & (abs(data_reco["paired_outer_eta"]) < 2.4) & (data_reco["Muon_pt_RC"] > 20) & (abs(data_reco["paired_inner_eta"]) < 2.5) 

     
    data_reco["criteria_probe_pass"] = data_reco["paired_Muon_ID"] >= 3
    data_reco["criteria_probe_fail"] = data_reco["paired_Muon_ID"] == 2
    
    
    # apply mask for all reconstructed muons:
    muon_plus = data_reco[["paired_outer_charge", "paired_outer_eta", "paired_outer_phi", "paired_outer_pt", "paired_nPV", "correct_charge_genPt", "paired_antiMuon_genEta", "paired_antiMuon_genPhi", "paired_antiMuon_genPt", "criteria_probe_pass", "criteria_probe_fail", "paired_criteria_tag"]][mask_plus]
    muon_minus = data_reco[["paired_outer_charge", "paired_outer_eta", "paired_outer_phi", "paired_outer_pt", "paired_nPV", "correct_charge_genPt", "paired_muon_genEta", "paired_muon_genPhi", "paired_muon_genPt", "criteria_probe_pass", "criteria_probe_fail", "paired_criteria_tag"]][mask_minus]


     
        
    # TO DO: apply corrections to outer, and its errors
        
    # scale_plus = rochester.kScaleDT(muon_plus["charge"], muon_plus["pt"], muon_plus["eta"], muon_plus["phi"])
    # scale_minus = rochester.kScaleDT(muon_minus["charge"], muon_minus["pt"], muon_minus["eta"], muon_minus["phi"])
    # spread_plus = rochester.kSpreadMC(muon_plus["charge"], muon_plus["pt"], muon_plus["eta"], muon_plus["phi"], muon_plus["correct_charge_genPt"])
    # spread_minus = rochester.kSpreadMC(muon_minus["charge"], muon_minus["pt"], muon_minus["eta"], muon_minus["phi"], muon_minus["correct_charge_genPt"])

    # scale_plus_error = rochester.kScaleDTerror(muon_plus["charge"], muon_plus["pt"], muon_plus["eta"], muon_plus["phi"])
    # scale_minus_error = rochester.kScaleDTerror(muon_minus["charge"], muon_minus["pt"], muon_minus["eta"], muon_minus["phi"])
    # spread_plus_error = rochester.kSpreadMCerror(muon_plus["charge"], muon_plus["pt"], muon_plus["eta"], muon_plus["phi"], muon_plus["correct_charge_genPt"])
    # spread_minus_error = rochester.kSpreadMCerror(muon_minus["charge"], muon_minus["pt"], muon_minus["eta"], muon_minus["phi"], muon_minus["correct_charge_genPt"])

   
    # all the templates
    # template_types = ["none", "nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]
    template_types = ["none"]
    result = {}


    muon_plus["Muon_pt_RC_none"] = muon_plus["paired_outer_pt"]
    muon_minus["Muon_pt_RC_none"] = muon_minus["paired_outer_pt"]
    
    # muon_plus["Muon_pt_RC_nominal"] = muon_plus["paired_outer_pt"] * scale_plus * spread_plus
    # muon_minus["Muon_pt_RC_nominal"] = muon_minus["paired_outer_pt"] * scale_minus * spread_minus

    # muon_plus["Muon_pt_RC_plus_scale"] = muon_plus["paired_outer_pt"] * (scale_plus + scale_plus_error) * spread_plus
    # muon_minus["Muon_pt_RC_plus_scale"] = muon_minus["paired_outer_pt"] * (scale_minus + scale_minus_error) * spread_minus

    # muon_plus["Muon_pt_RC_minus_scale"] = muon_plus["paired_outer_pt"] * (scale_plus - scale_plus_error) * spread_plus
    # muon_minus["Muon_pt_RC_minus_scale"] = muon_minus["paired_outer_pt"] * (scale_minus - scale_minus_error) * spread_minus

    # muon_plus["Muon_pt_RC_plus_spread"] = muon_plus["paired_outer_pt"] * scale_plus * (spread_plus + spread_plus_error)
    # muon_minus["Muon_pt_RC_plus_spread"] = muon_minus["paired_outer_pt"] * scale_minus * (spread_minus + spread_minus_error)

    # muon_plus["Muon_pt_RC_minus_spread"] = muon_plus["paired_outer_pt"] * scale_plus * (spread_plus - spread_plus_error)
    # muon_minus["Muon_pt_RC_minus_spread"] = muon_minus["paired_outer_pt"] * scale_minus * (spread_minus - spread_minus_error)


    for rochester_error_type in template_types:
            
        # mask for pt
        muon_sta_plus = muon_plus[muon_plus["Muon_pt_RC_{0}".format(rochester_error_type)] > 25]
        muon_sta_minus = muon_minus[muon_minus["Muon_pt_RC_{0}".format(rochester_error_type)] > 25]
            
        
        # define mask for Glopass:

        pairs_1_pass = ak.cartesian([muon_sta_plus["paired_criteria_tag"], muon_sta_minus["criteria_probe_pass"]])
        l_idx_1_pass, r_idx_1_pass = ak.unzip(pairs_1_pass)
        mask_1_pass = (l_idx_1_pass == r_idx_1_pass) & (l_idx_1_pass + r_idx_1_pass == True)

        pairs_2_pass = ak.cartesian([muon_sta_plus["criteria_probe_pass"], muon_sta_minus["paired_criteria_tag"]])
        l_idx_2_pass, r_idx_2_pass = ak.unzip(pairs_2_pass)
        mask_2_pass = (l_idx_2_pass == r_idx_2_pass) & (l_idx_2_pass + r_idx_2_pass == True)

        mask_pass = mask_1_pass + mask_2_pass

        # define mask for Glofail:

        pairs_1_fail = ak.cartesian([muon_sta_plus["paired_criteria_tag"], muon_sta_minus["criteria_probe_fail"]])
        l_idx_1_fail, r_idx_1_fail = ak.unzip(pairs_1_fail)
        mask_1_fail = (l_idx_1_fail == r_idx_1_fail) & (l_idx_1_fail + r_idx_1_fail == True)

        pairs_2_fail = ak.cartesian([muon_sta_plus["criteria_probe_fail"], muon_sta_minus["paired_criteria_tag"]])
        l_idx_2_fail, r_idx_2_fail = ak.unzip(pairs_2_fail)
        mask_2_fail = (l_idx_2_fail == r_idx_2_fail) & (l_idx_2_fail + r_idx_2_fail == True)

        mask_fail = mask_1_fail + mask_2_fail
               
            
        # pairing, and applying Glo masks:
        mask_types = ["Glopass", "Glofail"]
        for Glo_type in mask_types:
            if Glo_type == "Glopass":
                mask = mask_pass
            elif Glo_type == "Glofail":
                mask = mask_fail
            
            cuts = ["gen_cut", "reco_cut"]
            for cut in cuts:
                if cut == "gen_cut":
                    eta_values_plus, eta_values_minus, phi_values_plus, phi_values_minus, pt_values_plus, pt_values_minus = "paired_antiMuon_genEta", "paired_muon_genEta", "paired_antiMuon_genPhi", "paired_muon_genPhi", "paired_antiMuon_genPt", "paired_muon_genPt"
                elif cut == "reco_cut":
                    eta_values_plus, eta_values_minus, phi_values_plus, phi_values_minus, pt_values_plus, pt_values_minus = "paired_outer_eta", 'paired_outer_eta', "paired_outer_phi", "paired_outer_phi", "Muon_pt_RC_{0}".format(rochester_error_type), "Muon_pt_RC_{0}".format(rochester_error_type)


                # pairing, and applying masks:
                pairs_eta = ak.cartesian([muon_sta_plus[eta_values_plus], muon_sta_minus[eta_values_minus]])[mask]
                pairs_phi = ak.cartesian([muon_sta_plus[phi_values_plus], muon_sta_minus[phi_values_minus]])[mask]
                pairs_pt = ak.cartesian([muon_sta_plus[pt_values_plus], muon_sta_minus[pt_values_minus]])[mask]
                pairs_PV = ak.cartesian([muon_sta_plus["paired_nPV"], muon_sta_minus["paired_nPV"]])[mask]


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

                mass0 = (mu1+mu2).mass
                mask_repeats = np.concatenate((np.array([1]), np.diff(mass0))) != 0

                mass = mass0[mask_repeats]
                PV = l_PV[mask_repeats]


                # binning
                if bin_in_PV == True:
                    hist = Hist.new.Regular(70, 56, 126, name="mass").Regular(100, 0.5, 100.5, name="PV").Double()
                    hist.fill(mass, PV)
                    result["hist_{0}_{1}".format(rochester_error_type, Glo_type)] = hist.view()
                else:
                    hist = Hist.new.Regular(70, 56, 126, name="mass").Double()
                    hist.fill(mass)
                    result["hist_{0}_{1}".format(rochester_error_type, Glo_type)] =  hist.view()           
    
    return result


#####################################

def generate_templates_Sta(data_reco, rochester, bin_in_PV = True):

    # define outer pt, eta, phi, charge
    data_reco["outer_pt"] = data_reco["Muon_useUpdated"] * data_reco["Muon_ptStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_ptStaReg"] 
    data_reco["outer_eta"] = data_reco["Muon_useUpdated"] * data_reco["Muon_etaStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_etaStaReg"] 
    data_reco["outer_phi"] = data_reco["Muon_useUpdated"] * data_reco["Muon_phiStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_phiStaReg"] 
    data_reco["outer_charge"] = data_reco["Muon_useUpdated"] * data_reco["Muon_chargeStaUpd"] + (1 - data_reco["Muon_useUpdated"]) * data_reco["Muon_chargeStaReg"] 
    
       
    # define inner parameters
    inner_types1 = ["pt", "eta", "phi", "charge"]
    inner_types2 = ["nPixelHits", "nTrackerLayers", "trackAlgo"]
    for inner_type in inner_types1:
        data_reco["inner_{0}".format(inner_type)] = ak.concatenate((data_reco["Muon_{0}Trk".format(inner_type)], data_reco["Track_{0}".format(inner_type)]), axis = 1)
    for inner_type in inner_types2:
        data_reco["{0}".format(inner_type)] = ak.concatenate((data_reco["Muon_{0}".format(inner_type)], data_reco["Track_{0}".format(inner_type)]), axis = 1)

    
    data_reco2 = {}
    # pairing inner and outer tracks
    kinematics = ["pt", "eta", "phi", "charge"]
    for value in kinematics:
        data_reco2["pairing_{0}".format(value)] = ak.cartesian([data_reco["inner_{0}".format(value)], data_reco["outer_{0}".format(value)]])
        data_reco["paired_inner_{0}".format(value)], data_reco["paired_outer_{0}".format(value)] = ak.unzip(data_reco2["pairing_{0}".format(value)])

    # define delta_R2 between inner and outer
    data_reco["delR2_io"] = (data_reco["paired_outer_eta"] - data_reco["paired_inner_eta"])**2 + (data_reco["paired_outer_phi"] - data_reco["paired_inner_phi"])**2


    # define tag criteria
    data_reco["HLT_trigger_pass"] = data_reco["Muon_triggerBits"]&1 == 1
    data_reco["criteria_tag_muons"] = (data_reco["HLT_trigger_pass"] == True) & (data_reco["Muon_ID"] >= 4)
    data_reco["criteria_tag_tracks"] = ak.ones_like(data_reco["Track_pt"])==0
    data_reco["criteria_tag"] = ak.concatenate((data_reco["criteria_tag_muons"], data_reco["criteria_tag_tracks"]), axis = 1)

        
    # reshape gen-type variables
    variables_gen = ["nPV", "muon_genPt", "muon_genEta", "muon_genPhi", "antiMuon_genPt", "antiMuon_genEta", "antiMuon_genPhi"]
    for variable in variables_gen:
        data_reco["paired_{0}".format(variable)]  = ak.ones_like(data_reco["delR2_io"]) * data_reco["{0}".format(variable)]     
        
    # reshape nongen-type variables
    variables_nongen = ["nPixelHits", "nTrackerLayers", "trackAlgo", "criteria_tag"]
    for variable in variables_nongen:
        data_reco2["pairing_{0}".format(variable)] = ak.cartesian([data_reco["{0}".format(variable)], data_reco["outer_pt".format(variable)]])
        data_reco["paired_{0}".format(variable)], data_reco2["paired_{0}".format(variable)] = ak.unzip(data_reco2["pairing_{0}".format(variable)])

    
    
    # define delta_R2 with gen
    data_reco["delR2_minus"] = (data_reco["paired_inner_eta"] - data_reco["muon_genEta"])**2 + (data_reco["paired_inner_phi"] - data_reco["muon_genPhi"])**2
    data_reco["delR2_plus"] = (data_reco["paired_inner_eta"] - data_reco["antiMuon_genEta"])**2 + (data_reco["paired_inner_phi"] - data_reco["antiMuon_genPhi"])**2
        
        
    data_reco["correct_charge_genPt"] = (data_reco["paired_inner_charge"] + 1)/2 * data_reco["antiMuon_genPt"] - (data_reco["paired_inner_charge"] - 1)/2 * data_reco["muon_genPt"]

    # define masks
    mask_minus = (data_reco["delR2_minus"]<0.09) & (data_reco["paired_inner_charge"] == -1) & (abs(data_reco["paired_inner_eta"]) < 2.4) & (data_reco["paired_inner_pt"] > 25) & (data_reco["paired_nPixelHits"] > 0 ) & (data_reco["paired_nTrackerLayers"] > 5) & (data_reco["paired_trackAlgo"] != 13) & (data_reco["paired_trackAlgo"] != 14) & (data_reco["paired_outer_pt"] > 20) & (abs(data_reco["paired_outer_eta"]) < 2.5) 
    mask_plus  = (data_reco["delR2_plus"]<0.09) & (data_reco["paired_inner_charge"] == 1)  & (abs(data_reco["paired_inner_eta"]) < 2.4) & (data_reco["paired_inner_pt"] > 25) & (data_reco["paired_nPixelHits"] > 0 ) & (data_reco["paired_nTrackerLayers"] > 5) & (data_reco["paired_trackAlgo"] != 13) & (data_reco["paired_trackAlgo"] != 14) & (data_reco["paired_outer_pt"] > 20) & (abs(data_reco["paired_outer_eta"]) < 2.5) 

        
    data_reco["criteria_probe_pass"] = (data_reco["delR2_io"]<0.09)
    data_reco["criteria_probe_fail"] = (data_reco["delR2_io"]>=0.09)
        
        
    muon_plus = data_reco[["paired_inner_charge", "paired_inner_eta", "paired_inner_phi", "paired_inner_pt", "paired_nPV", "correct_charge_genPt", "paired_antiMuon_genEta", "paired_antiMuon_genPhi", "paired_antiMuon_genPt", "criteria_probe_pass", "criteria_probe_fail", "paired_criteria_tag"]][mask_plus]
    muon_minus = data_reco[["paired_inner_charge", "paired_inner_eta", "paired_inner_phi", "paired_inner_pt", "paired_nPV", "correct_charge_genPt", "paired_muon_genEta", "paired_muon_genPhi", "paired_muon_genPt", "criteria_probe_pass", "criteria_probe_fail", "paired_criteria_tag"]][mask_minus]


    # rochester corrections
    scale_plus = rochester.kScaleDT(muon_plus["paired_inner_charge"], muon_plus["paired_inner_pt"], muon_plus["paired_inner_eta"], muon_plus["paired_inner_phi"])
    scale_minus = rochester.kScaleDT(muon_minus["paired_inner_charge"], muon_minus["paired_inner_pt"], muon_minus["paired_inner_eta"], muon_minus["paired_inner_phi"])
    spread_plus = rochester.kSpreadMC(muon_plus["paired_inner_charge"], muon_plus["paired_inner_pt"], muon_plus["paired_inner_eta"], muon_plus["paired_inner_phi"], muon_plus["correct_charge_genPt"])
    spread_minus = rochester.kSpreadMC(muon_minus["paired_inner_charge"], muon_minus["paired_inner_pt"], muon_minus["paired_inner_eta"], muon_minus["paired_inner_phi"], muon_minus["correct_charge_genPt"])

    scale_plus_error = rochester.kScaleDTerror(muon_plus["paired_inner_charge"], muon_plus["paired_inner_pt"], muon_plus["paired_inner_eta"], muon_plus["paired_inner_phi"])
    scale_minus_error = rochester.kScaleDTerror(muon_minus["paired_inner_charge"], muon_minus["paired_inner_pt"], muon_minus["paired_inner_eta"], muon_minus["paired_inner_phi"])
    spread_plus_error = rochester.kSpreadMCerror(muon_plus["paired_inner_charge"], muon_plus["paired_inner_pt"], muon_plus["paired_inner_eta"], muon_plus["paired_inner_phi"], muon_plus["correct_charge_genPt"])
    spread_minus_error = rochester.kSpreadMCerror(muon_minus["paired_inner_charge"], muon_minus["paired_inner_pt"], muon_minus["paired_inner_eta"], muon_minus["paired_inner_phi"], muon_minus["correct_charge_genPt"])

        
    # all the rochester templates
    template_types = ["none", "nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]

    muon_plus["Muon_pt_RC_none"] = muon_plus["paired_inner_pt"]
    muon_minus["Muon_pt_RC_none"] = muon_minus["paired_inner_pt"]

    muon_plus["Muon_pt_RC_nominal"] = muon_plus["paired_inner_pt"] * scale_plus * spread_plus
    muon_minus["Muon_pt_RC_nominal"] = muon_minus["paired_inner_pt"] * scale_minus * spread_minus

    muon_plus["Muon_pt_RC_plus_scale"] = muon_plus["paired_inner_pt"] * (scale_plus + scale_plus_error) * spread_plus
    muon_minus["Muon_pt_RC_plus_scale"] = muon_minus["paired_inner_pt"] * (scale_minus + scale_minus_error) * spread_minus

    muon_plus["Muon_pt_RC_minus_scale"] = muon_plus["paired_inner_pt"] * (scale_plus - scale_plus_error) * spread_plus
    muon_minus["Muon_pt_RC_minus_scale"] = muon_minus["paired_inner_pt"] * (scale_minus - scale_minus_error) * spread_minus

    muon_plus["Muon_pt_RC_plus_spread"] = muon_plus["paired_inner_pt"] * scale_plus * (spread_plus + spread_plus_error)
    muon_minus["Muon_pt_RC_plus_spread"] = muon_minus["paired_inner_pt"] * scale_minus * (spread_minus + spread_minus_error)

    muon_plus["Muon_pt_RC_minus_spread"] = muon_plus["paired_inner_pt"] * scale_plus * (spread_plus - spread_plus_error)
    muon_minus["Muon_pt_RC_minus_spread"] = muon_minus["paired_inner_pt"] * scale_minus * (spread_minus - spread_minus_error)
    
    
    result = {}
    for rochester_error_type in template_types:

        # mask for pt
        muon_trk_plus = muon_plus[muon_plus["Muon_pt_RC_{0}".format(rochester_error_type)] > 25]
        muon_trk_minus = muon_minus[muon_minus["Muon_pt_RC_{0}".format(rochester_error_type)] > 25]


        # define mask for Stapass:
        pairs_1_pass = ak.cartesian([muon_trk_plus["paired_criteria_tag"], muon_trk_minus["criteria_probe_pass"]])
        l_idx_1_pass, r_idx_1_pass = ak.unzip(pairs_1_pass)
        mask_1_pass = (l_idx_1_pass == r_idx_1_pass) & (l_idx_1_pass + r_idx_1_pass == True)

        pairs_2_pass = ak.cartesian([muon_trk_plus["criteria_probe_pass"], muon_trk_minus["paired_criteria_tag"]])
        l_idx_2_pass, r_idx_2_pass = ak.unzip(pairs_2_pass)
        mask_2_pass = (l_idx_2_pass == r_idx_2_pass) & (l_idx_2_pass + r_idx_2_pass == True)

        mask_pass = mask_1_pass + mask_2_pass

        # define mask for Stafail:
        pairs_1_fail = ak.cartesian([muon_trk_plus["paired_criteria_tag"], muon_trk_minus["criteria_probe_fail"]])
        l_idx_1_fail, r_idx_1_fail = ak.unzip(pairs_1_fail)
        mask_1_fail = (l_idx_1_fail == r_idx_1_fail) & (l_idx_1_fail + r_idx_1_fail == True)

        pairs_2_fail = ak.cartesian([muon_trk_plus["criteria_probe_fail"], muon_trk_minus["paired_criteria_tag"]])
        l_idx_2_fail, r_idx_2_fail = ak.unzip(pairs_2_fail)
        mask_2_fail = (l_idx_2_fail == r_idx_2_fail) & (l_idx_2_fail + r_idx_2_fail == True)

        mask_fail = mask_1_fail + mask_2_fail


        # pairing, and applying Sta masks:
        mask_types = ["Stapass", "Stafail"]
        for Sta_type in mask_types:
            if Sta_type == "Stapass":
                mask = mask_pass 
            elif Sta_type == "Stafail":
                mask = mask_fail


            cuts = ["gen_cut", "reco_cut"]
            for cut in cuts:

                if cut == "gen_cut":
                    eta_values_plus, eta_values_minus, phi_values_plus, phi_values_minus, pt_values_plus, pt_values_minus = "paired_antiMuon_genEta", "paired_muon_genEta", "paired_antiMuon_genPhi", "paired_muon_genPhi", "paired_antiMuon_genPt", "paired_muon_genPt"
                elif cut == "reco_cut":
                    eta_values_plus, eta_values_minus, phi_values_plus, phi_values_minus, pt_values_plus, pt_values_minus = "paired_inner_eta", 'paired_inner_eta', "paired_inner_phi", "paired_inner_phi", "Muon_pt_RC_{0}".format(rochester_error_type), "Muon_pt_RC_{0}".format(rochester_error_type)


                # pairing, and applying masks:
                pairs_eta = ak.cartesian([muon_trk_plus[eta_values_plus], muon_trk_minus[eta_values_minus]])[mask]
                pairs_phi = ak.cartesian([muon_trk_plus[phi_values_plus], muon_trk_minus[phi_values_minus]])[mask]
                pairs_pt = ak.cartesian([muon_trk_plus[pt_values_plus], muon_trk_minus[pt_values_minus]])[mask]
                pairs_PV = ak.cartesian([muon_trk_plus["paired_nPV"], muon_trk_minus["paired_nPV"]])[mask]

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

                mass0 = (mu1+mu2).mass
                mask_repeats = np.concatenate((np.array([1]), np.diff(mass0))) != 0

                mass = mass0[mask_repeats]
                PV = l_PV[mask_repeats]


                # binning
                if bin_in_PV == True:
                    hist = Hist.new.Regular(50, 66, 116, name="mass").Regular(100, 0.5, 100.5, name="PV").Double()
                    hist.fill(mass, PV)
                    result["hist_{0}_{1}_{2}".format(rochester_error_type, Sta_type, cut)] = hist.view()
                else:
                    hist = Hist.new.Regular(50, 66, 116, name="mass").Double()
                    hist.fill(mass)
                    result["hist_{0}_{1}_{2}".format(rochester_error_type, Sta_type, cut)] =  hist.view()

                    
    return result

result = generate_templates_Sta(data_reco, rochester)


#####################################



parser = argparse.ArgumentParser()
parser.add_argument("mc_data_file")
parser.add_argument("rochester_file")
args = parser.parse_args()


filename_ntuples = args.mc_data_file
rochester_filename = args.rochester_file


# filename_ntuples = "C:/Users/byron/ZCounting/DYJetsToLL_M_50_LO_FlatPU0to75_Autumn18.root"
# rochester_filename = "C:/Users/byron/ZCounting/ZHarvester/res/Rocco/RoccoR2018.txt"


f1 = uproot.open(filename_ntuples)

treename = "zcounting/tree"
data = uproot.open(filename_ntuples+":"+treename)

data_reco = data.arrays(["muon_genPt","muon_genEta", "muon_genPhi",
 "antiMuon_genPt","antiMuon_genEta","antiMuon_genPhi",
 "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge",
 "Muon_ID", "Muon_triggerBits", "nPV",
 "Muon_useUpdated", "Muon_ptStaUpd", "Muon_etaStaUpd", "Muon_phiStaUpd", "Muon_chargeStaUpd",
 "Muon_ptStaReg", "Muon_etaStaReg", "Muon_phiStaReg", "Muon_chargeStaReg",
 "Track_pt", "Track_eta", "Track_phi", "Track_charge",
 "Muon_ptTrk", "Muon_etaTrk", "Muon_phiTrk", "Muon_chargeTrk",
 "Track_nPixelHits", "Track_nTrackerLayers", "Track_trackAlgo",
 "Muon_nPixelHits", "Muon_nTrackerLayers","Muon_trackAlgo"], "(decayMode==13) & (nMuon>=2)")

# data_reco = data.arrays(["muon_genPt","muon_genEta", "muon_genPhi",
#  "antiMuon_genPt","antiMuon_genEta","antiMuon_genPhi",
#  "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge",
#  "Muon_ID", "Muon_triggerBits", "nPV"], "(decayMode==13) & (nMuon>=2)")


rochester = _get_rochester_corr(rochester_filename)

templates = generate_templates_Sta(data_reco, rochester, bin_in_PV = True)[0]
# rochester_corrections = generate_templates_Glo(data_reco, rochester)[1]


with h5py.File('templates_Sta.hdf5', 'w') as outfile:
    for dset_name in templates:
        dset = outfile.create_dataset(dset_name, data = templates[dset_name])

# with h5py.File('rochester_IDfail.hdf5', 'w') as outfile:
#     for dset_name in rochester_corrections:
#         dset = outfile.create_dataset(dset_name, data = rochester_corrections[dset_name])


# pdb.set_trace()