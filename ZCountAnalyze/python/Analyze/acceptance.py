#!/usr/bin/env python3
#
################################################################################
#
# calculate the acceptance factor for a given selection on gen level objects
#
################################################################################

import uproot 
import pandas as pd
import numpy as np
import pdb
import uncertainties as unc
import awkward as ak
import vector

### acceptance cuts
ptCut = 25
etaCut = 2.4
mass_lo = 60
mass_hi = 120
final_states = "electrons"

if final_states == "muons":
    decay_mode = 13
    mass = 0.105658369
elif final_states == "electrons":
    decay_mode = 11
    mass = 0.000510998

undressed = True
dressed = True

### resources
filename_ntuples = "/eos/cms/store/group/comm_luminosity/ZCounting/2022/MC/output/V01/DYJetsToLL_M_50_LO_FlatPU0to70_Run3Winter22.root"
treename = "zcounting/tree"

branches = ["lepton_genPt","lepton_genEta","lepton_genPhi",
    "antiLepton_genPt","antiLepton_genEta","antiLepton_genPhi", 
    "z_genMass", "eventweight", "MEWeight", "PSWeight"]

branches_dressed = [ 
    "GenDressedLepton_pt", "GenDressedLepton_eta", "GenDressedLepton_phi", "GenDressedLepton_pdgId",
    ]
branches_weights = ["eventweight", "MEWeight", "PSWeight"]

aliases = {

}

### eventweights to estimate uncertainty 
# TODO: check if the indices for ME weights are correct
weights = {
    "muRmuF_up":   "MEWeight[0]",
    "muRmuF_down": "MEWeight[1]",
    "muR_up":      "MEWeight[2]",
    "muR_down":    "MEWeight[3]",
    "muF_up":      "MEWeight[4]",
    "muF_down":    "MEWeight[5]",
    "ISR_up":      "PSWeight[0]",
    "ISR_down":    "PSWeight[2]",
    "FSR_up":      "PSWeight[1]",
    "FSR_down":    "PSWeight[3]",
}

n_total = {"nominal" : unc.ufloat(0.,0.)}
n_selected = {"nominal" : unc.ufloat(0.,0.)}

for key in weights.keys():
    n_total[key] = 0
    n_selected[key] = 0

print(f"open tree `{treename}` from file `{filename_ntuples}`")
# with uproot.open(filename_ntuples+":"+treename) as events:
    
if undressed:
    for batch in uproot.iterate(filename_ntuples+":"+treename, step_size="100 MB", library="ak", 
        filter_name=branches, aliases=aliases,
        cut=f"(decayMode=={decay_mode})"# & (z_genMass > {mass_lo}) & (z_genMass < {mass_hi})"
    ):
        
        print(f"Have {len(batch)} events")
        
        # acceptance cuts
        print(f"Select events in acceptance")
        selected = batch[
            (batch["lepton_genPt"] > ptCut) & (batch["antiLepton_genPt"] > ptCut) 
            & (abs(batch["lepton_genEta"]) < etaCut) & (abs(batch["antiLepton_genEta"]) < etaCut)
            & (batch["z_genMass"] < mass_hi) & (batch["z_genMass"] > mass_lo)
            ]

        print(f"Compute expected event numbers")            
        n_total["nominal"] += unc.ufloat(sum(batch["eventweight"]), np.sqrt(sum((batch["eventweight"])**2)))
        n_selected["nominal"] += unc.ufloat(sum(selected["eventweight"]), np.sqrt(sum((selected["eventweight"])**2)))
            
        for key, weight in weights.items():
            if len(weight.split("[")) == 2:
                weight, index = weight.split("[")
                index = int(index.split("]")[0])
                n_total[key] += sum(batch[weight][:,index]*batch["eventweight"])
                n_selected[key] += sum(selected[weight][:,index]*selected["eventweight"])
            else:
                n_total[key] += sum(batch[weight]*batch["eventweight"])
                n_selected[key] += sum(selected[weight]*selected["eventweight"])
                
        print(f"Read next batch of events")

    # nominal acceptance
    acceptance = n_selected["nominal"] / n_total["nominal"]

    print("###################################################################")
    print(f"# Calculate acceptance for {final_states} with:")
    print(f"# pT > {ptCut}") 
    print(f"# |eta| < {etaCut}")
    print(f"# {mass_lo} < |mll| < {mass_hi} GeV")
    print(f"# A = {acceptance}")
    acceptances = {}
    for key in weights.keys():    
        acceptances[key] = n_selected[key] / n_total[key]
        print(f"# A({key}) = {acceptances[key]}")
    print("###################################################################")


if dressed:
    print("For dressed leptons:")

    for batch in uproot.iterate(filename_ntuples+":"+treename, step_size="100 MB", library="ak", 
        filter_name=branches_dressed+branches_weights, aliases=aliases,
        cut=f"(decayMode=={decay_mode})"# & (z_genMass > {mass_lo}) & (z_genMass < {mass_hi})"
    ):
        
        print(f"Have {len(batch)} events")
        
        # acceptance cuts
        print(f"Select events in acceptance")

        leptons = batch[branches_dressed]
        leptons = leptons[
            (leptons["GenDressedLepton_pt"] > ptCut)
            & (abs(leptons["GenDressedLepton_eta"]) < etaCut) 
            & (abs(leptons["GenDressedLepton_pdgId"]) == decay_mode)
            ]
        
        # require pairs with opposite charge
        p_charge = ak.combinations(leptons["GenDressedLepton_pdgId"], 2)        
        lefts, rights = ak.unzip(p_charge)
        mask_charge = lefts/decay_mode*rights/decay_mode == -1

        p_pt = ak.combinations(leptons["GenDressedLepton_pt"], 2)
        p_eta = ak.combinations(leptons["GenDressedLepton_eta"], 2)
        p_phi = ak.combinations(leptons["GenDressedLepton_phi"], 2)

        # 2.) calculate mass of each pair
        l_pt, r_pt = ak.unzip(p_pt)
        l_eta, r_eta = ak.unzip(p_eta)
        l_phi, r_phi = ak.unzip(p_phi)

        mu1 = vector.obj(pt=l_pt, phi=l_phi, eta=l_eta, mass=l_pt*0+mass)
        mu2 = vector.obj(pt=r_pt, phi=r_phi, eta=r_eta, mass=r_pt*0+mass)

        masses = (mu1 + mu2).mass

        mask_mass = (masses > mass_lo) & (masses < mass_hi)

        events = batch[branches_weights]
        events["mass"] = masses[mask_mass & mask_charge]

        selected = events[(ak.num(leptons["GenDressedLepton_pt"]) >=2) & (ak.num(events["mass"]) >=1) ]

        print(f"Compute expected event numbers")            
        n_total["nominal"] += unc.ufloat(sum(events["eventweight"]), np.sqrt(sum((events["eventweight"])**2)))
        n_selected["nominal"] += unc.ufloat(sum(selected["eventweight"]), np.sqrt(sum((selected["eventweight"])**2)))
            
        for key, weight in weights.items():
            if len(weight.split("[")) == 2:
                weight, index = weight.split("[")
                index = int(index.split("]")[0])
                n_total[key] += sum(events[weight][:,index]*events["eventweight"])
                n_selected[key] += sum(selected[weight][:,index]*selected["eventweight"])
            else:
                n_total[key] += sum(events[weight]*events["eventweight"])
                n_selected[key] += sum(selected[weight]*selected["eventweight"])
                
        print(f"Read next batch of events")

    # nominal acceptance
    acceptance = n_selected["nominal"] / n_total["nominal"]

    print("###################################################################")
    print(f"# Calculate acceptance for {final_states} with:")
    print(f"# pT > {ptCut}") 
    print(f"# |eta| < {etaCut}")
    print(f"# {mass_lo} < |mll| < {mass_hi} GeV")
    print(f"# A = {acceptance}")
    acceptances = {}
    for key in weights.keys():    
        acceptances[key] = n_selected[key] / n_total[key]
        print(f"# A({key}) = {acceptances[key]}")
    print("###################################################################")