import awkward as ak
import numpy as np

import h5py

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("templates_gen") # hdf5 file
parser.add_argument("templates_reco") # hdf5 file
args = parser.parse_args()


# load gen-cut histograms
templates_gen_hdf5 = h5py.File(args.templates_gen, 'r')

templates_gen = {}
for key in templates_gen_hdf5.keys():
    templates_gen[key] = np.array(templates_gen_hdf5[key])

# reduce to 1d template in PV
N_gen = {}
for key in templates_gen.keys():
    N_gen[key] = np.sum(templates_gen[key], axis = 0)


# just some useful quantities for later
alpha1 = 2 * N_gen["HLT2"] + N_gen["HLT1"]
alpha2 = alpha1 + N_gen["IDfail"]
alpha3 = alpha2 + N_gen["Glofail"]
beta3 = N_gen["HLT_Trkpass"] + N_gen["HLT_Trkfail"] + N_gen["ID_Trkfail"] + N_gen["Glo_Trkfail"]
beta4 = beta3 + N_gen["Sta_Trkfail"]


# calculate gen efficiencies
eff = {}
eff["HLT"] = (2 * N_gen["HLT2"]) / alpha1
eff["ID"] = alpha1 / alpha2
eff["Glo"] = N_gen["Glopass"] / (N_gen["Glopass"] + N_gen["Glofail"])
eff["Sta"] = N_gen["Sta_Trkpass"] / (N_gen["Sta_Trkpass"] + N_gen["Sta_Trkfail"])

# calculate correlations
correlations = {}

correlations["HLT"] = 4 * (N_gen["HLT2"] + N_gen["HLT1"] + N_gen["HLT0"]) * N_gen["HLT2"] / alpha1**2
correlations["ID"] = (N_gen["ID2"] + N_gen["ID1"] + N_gen["ID0"]) / (N_gen["HLT2"] + N_gen["HLT1"] + N_gen["HLT0"])  *  alpha1**2 / alpha2**2
correlations["Glo"] = N_gen["all"] / (N_gen["ID2"] + N_gen["ID1"] + N_gen["ID0"]) * alpha2**2 / alpha3**2 * beta3**2 / beta4**2
correlations["ID_total"] = correlations["ID"] * correlations["Glo"]


correlations["io"] = alpha2 / eff["Glo"] / eff["Sta"] / beta4 * beta3 / alpha3


# load reco-cut histograms
templates_reco_hdf5 = h5py.File(args.templates_reco, 'r')

templates_reco = {}
for key in templates_reco_hdf5.keys():
    templates_reco[key] = np.array(templates_reco_hdf5[key])

# reduce to 1d template in PV
N_reco = {}
for key in templates_reco.keys():
    N_reco[key] = np.sum(templates_reco[key], axis = 0)


N = {}
for template_type in ["HLT2", "HLT1", "IDfail", "Stapass", "Stafail"]: # need to include for Glopass and Glofail when we have corrections for them
    N["{0}_nominal".format(template_type)] = N_reco["{0}_nominal".format(template_type)]
    N["{0}_plusalpha".format(template_type)] = N_reco["{0}_plus_scale".format(template_type)]
    N["{0}_minusalpha".format(template_type)] = N_reco["{0}_minus_scale".format(template_type)]
    N["{0}_plusbeta".format(template_type)] = N_reco["{0}_plus_spread".format(template_type)]
    N["{0}_minusbeta".format(template_type)] = N_reco["{0}_minus_spread".format(template_type)]


acceptance = {}
acceptance["nominal"] = N_gen["all"] / (2 * N["HLT2_nominal"] + N["HLT1_nominal"] + N["IDfail_nominal"])**2 * (4 * N["HLT2_nominal"] * (N_reco["Glopass_none"] / (N_reco["Glopass_none"] + N_reco["Glofail_none"]) * N["Stapass_nominal"] / (N["Stapass_nominal"] + N["Stafail_nominal"]))**2) \
 / correlations["HLT"] / correlations["ID_total"] * correlations["io"]**2

correction_values = ["plusalpha", "minusalpha", "plusbeta", "minusbeta"]
for value in correction_values:
    acceptance["HLT2_{0}".format(value)] = N_gen["all"] / (2 * N["HLT2_{0}".format(value)] + N["HLT1_nominal"] + N["IDfail_nominal"])**2 * (4 * N["HLT2_{0}".format(value)] * (N_reco["Glopass_none"] / (N_reco["Glopass_none"] + N_reco["Glofail_none"]) * N["Stapass_nominal"] / (N["Stapass_nominal"] + N["Stafail_nominal"]))**2) \
     / correlations["HLT"] / correlations["ID_total"] * correlations["io"]**2

for value in correction_values:
    acceptance["HLT1_{0}".format(value)] = acceptance["nominal"] = N_gen["all"] / (2 * N["HLT2_nominal"] + N["HLT1_{0}".format(value)] + N["IDfail_nominal"])**2 * (4 * N["HLT2_nominal"] * (N_reco["Glopass_none"] / (N_reco["Glopass_none"] + N_reco["Glofail_none"]) * N["Stapass_nominal"] / (N["Stapass_nominal"] + N["Stafail_nominal"]))**2) \
     / correlations["HLT"] / correlations["ID_total"] * correlations["io"]**2

for value in correction_values:
    acceptance["IDfail_{0}".format(value)] = N_gen["all"] / (2 * N["HLT2_nominal"] + N["HLT1_nominal"] + N["IDfail_{0}".format(value)])**2 * (4 * N["HLT2_nominal"] * (N_reco["Glopass_none"] / (N_reco["Glopass_none"] + N_reco["Glofail_none"]) * N["Stapass_nominal"] / (N["Stapass_nominal"] + N["Stafail_nominal"]))**2) \
     / correlations["HLT"] / correlations["ID_total"] * correlations["io"]**2

#include here for Glopass and Glofail below when ready

for value in correction_values:
    acceptance["Stapass_{0}".format(value)] = N_gen["all"] / (2 * N["HLT2_nominal"] + N["HLT1_nominal"] + N["IDfail_nominal"])**2 * (4 * N["HLT2_nominal"] * (N_reco["Glopass_none"] / (N_reco["Glopass_none"] + N_reco["Glofail_none"]) * N["Stapass_{0}".format(value)] / (N["Stapass_{0}".format(value)] + N["Stafail_nominal"]))**2) \
     / correlations["HLT"] / correlations["ID_total"] * correlations["io"]**2

for value in correction_values:
    acceptance["Stafail_{0}".format(value)] = N_gen["all"] / (2 * N["HLT2_nominal"] + N["HLT1_nominal"] + N["IDfail_nominal"])**2 * (4 * N["HLT2_nominal"] * (N_reco["Glopass_none"] / (N_reco["Glopass_none"] + N_reco["Glofail_none"]) * N["Stapass_nominal"] / (N["Stapass_nominal"] + N["Stafail_{0}".format(value)]))**2) \
     / correlations["HLT"] / correlations["ID_total"] * correlations["io"]**2


with h5py.File('templates/correlations.hdf5', 'w') as outfile:
    for dset_name in correlations:
        dset = outfile.create_dataset(dset_name, data = correlations[dset_name])

with h5py.File('templates/acceptance.hdf5', 'w') as outfile:
    for dset_name in acceptance:
        dset = outfile.create_dataset(dset_name, data = acceptance[dset_name])
