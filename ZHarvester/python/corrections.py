import json
import pdb
import os
import pickle

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import exp, pexp, plinear, plinear_step

def apply_muon_prefire(df):

    print("apply muon prefire corrections")

    # --- get prefire corrections
    with open('res/prefiring.json') as file_prefire:
        res_prefire = json.load(file_prefire)

    for lo, hi, era in (
        (272007, 278769, "2016preVFP"),
        (278769, 280919, "2016postVFP"),
        (280919, 284045, "2016H"),
        (297020, 307083, "2017"),
        (315252, 325274, "2018"),
        ):

        loc = (df['run'] >= lo) & (df['run'] < hi)

        for region in ("BB","BE","EE"):
            factor = float(res_prefire[era]["prefMuon"][region]['nominal'])
            df.loc[loc,'zDel{0}_mc'.format(region)] *= factor

def apply_pileup_correction(df):
    print("apply pileup corrections")

    # --- get pileup corrections

    for lo, hi, era in (
        (272007, 278769, "2016preVFP"),
        (278769, 284045, "2016postVFP"),
        (297020, 307083, "2017"),
        (315252, 325274, "2018"),
        ):

        loc = (df['run'] >= lo) & (df['run'] < hi)

        if sum(loc) == 0:
            continue

        # --- get pileup corrections
        with open('res/mcCorrections/corrections_nPU_{0}.p'.format(era), "r") as file_pileup:
            functions = pickle.load(file_pileup)

        for region in ("BB","BE","EE"):
            function = functions["effReco"+region]
            params = function['params']
            name = function['name']
            if name == "linear_step":
                df.loc[loc, 'zDel{0}_mc'.format(region)] *= plinear_step(df.loc[loc,'pileUp'], *params)
            elif name == "linear":
                df.loc[loc, 'zDel{0}_mc'.format(region)] *= plinear(df.loc[loc,'pileUp'], *params)


def apply_ECAL_prefire(df):
    print("apply ECAL prefire corrections")

    # --- get prefire corrections
    with open('res/prefiring.json') as file_prefire:
        res_prefire = json.load(file_prefire)

    for lo, hi, era in (
        (272007, 278769, "2016preVFP"),
        (278769, 284045, "2016postVFP"),
        (297020, 307083, "2017")
        ):

        loc = (df['run'] >= lo) & (df['run'] < hi)

        for region in ("BB","BE","EE"):
            params = res_prefire[era]["prefECAL"][region]['nominal']
            df.loc[loc, 'zDel{0}_mc'.format(region)] *= pexp(df.loc[loc,'pileUp'], *params)
