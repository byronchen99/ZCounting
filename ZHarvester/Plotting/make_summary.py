import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
import uncertainties as unc
import pdb
import matplotlib.pyplot as plt
import json

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime, unorm
from ZUtils.python.utils import pexp, pquad

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("-r","--rates", required=True, nargs='+', help="Nominator csv file with z rates per lumisection")
parser.add_argument("-x","--xsec",  required=True, type=str,
    help="csv file with z rates per lumisection where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Data Acquisition ##########

# --- get prefire corrections
prefire_variations_Muon = ('nominal', 'SystUp', 'SystDown', 'StatUp', 'StatDown')
prefire_variations_ECAL = ("nominal", "Up", "Down",)
with open('res/prefiring.json') as file_prefire:
    res_prefire = json.load(file_prefire)

# --- get Z xsec
print("get Z cross section")
data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)

# --- z luminosity
print("get Z luminosity")
data = pd.concat([pd.read_csv(csv, sep=',',low_memory=False) for csv in args.rates] +[data_xsec,], ignore_index=True, sort=False)

data['zRec_mc'] = data['zRecBB_mc'] + data['zRecBE_mc'] + data['zRecEE_mc']
data['zDel_mc'] = data['zDelBB_mc'] + data['zDelBE_mc'] + data['zDelEE_mc']
data['zDel_mcUp'] = data['zDelBB_mcUp'] + data['zDelBE_mcUp'] + data['zDelEE_mcUp']
# data['zDel_mcDown'] = data['zDelBB_mcDown'] + data['zDelBE_mcDown'] + data['zDelEE_mcDown']

# --->>> prefire corrections
print("apply muon prefire corrections")
for var in prefire_variations_Muon:
    for region in ("BB","BE","EE","I"):

        data['zDel{0}_mc_prefMuon_{1}'.format(region,var)] = data['zDel{0}_mc'.format(region)]

        for lo, hi, era in (
            (272007, 278769, "2016preVFP"),
            (278769, 280919, "2016postVFP"),
            (280919, 284045, "2016H"),
            (297020, 306463, "2017"),
            (306828, 307083, "2017"),   # 2017H
            (315252, 325274, "2018"),
        ):
            factor = float(res_prefire[era]["prefMuon"][region][var])
            loc = (data['run'] >= lo) & (data['run'] < hi)
            data.loc[loc,'zDel{0}_mc_prefMuon_{1}'.format(region,var)] *= factor

    data['zDel_mc_prefMuon_'+var] = data['zDelBB_mc_prefMuon_'+var] \
        + data['zDelBE_mc_prefMuon_'+var] \
        + data['zDelEE_mc_prefMuon_'+var]

print("apply ECAL prefire corrections")

for var in prefire_variations_ECAL:
    for region in ("BB","BE","EE","I"):

        data['zDel{0}_mc_prefECAL_{1}'.format(region,var)] = data['zDel{0}_mc_prefMuon_nominal'.format(region)]

        for lo, hi, era, func in (
            (272007, 278769, "2016preVFP", pexp),
            (278769, 284045, "2016postVFP", pexp),
            (297020, 306463, "2017", pexp),
            (306828, 307083, "2017H", pquad),   # 2017H
        ):
            params = res_prefire[era]["prefECAL"][region][var]
            loc = (data['run'] >= lo) & (data['run'] < hi)
            data.loc[loc, 'zDel{0}_mc_prefECAL_{1}'.format(region,var)] *= func(data.loc[loc,'pileUp'], *params)

    data['zDel_mc_prefECAL_'+var] = data['zDelBB_mc_prefECAL_'+var] \
        + data['zDelBE_mc_prefECAL_'+var] \
        + data['zDelEE_mc_prefECAL_'+var]

print("apply prefire corrections - done")
# <<<---

for region in ("","BB","BE","EE","I"):
    info = {}

    for lo, hi, era in (
        (272007, 278769, "2016preVFP"),
        (278769, 280919, "2016postVFP"),
        (280919, 284045, "2016H"),
        (297020, 306463, "2017"),
        (306828, 307083, "2017H"),   # 2017H
        (315252, 325274, "2018"),
    ):
        info[era] = {}
        info[era]['zRec_mc'] = sum(data.loc[(data['run'] >= lo) & (data['run'] < hi),'zRec{0}_mc'.format(region)])
        # info[era]['raw'] = sum(data.loc[(data['run'] >= lo) & (data['run'] < hi),'zDel'])
        info[era]['mc'] = sum(data.loc[(data['run'] >= lo) & (data['run'] < hi),'zDel{0}_mc'.format(region)])
        info[era]['mcUp'] = sum(data.loc[(data['run'] >= lo) & (data['run'] < hi),'zDel{0}_mcUp'.format(region)])
        # info[era]['mcDown'] = sum(data.loc[(data['run'] >= lo) & (data['run'] < hi),'zDel{0}_mcDown'.format(region)])
        for var in prefire_variations_Muon:
            info[era]["prefMuon_"+var] = sum(data.loc[(data['run'] >= lo) & (data['run'] < hi),'zDel{0}_mc_prefMuon_{1}'.format(region, var)])
        for var in prefire_variations_ECAL:
            info[era]["prefECAL_"+var] = sum(data.loc[(data['run'] >= lo) & (data['run'] < hi),'zDel{0}_mc_prefECAL_{1}'.format(region, var)])

    with open(outDir+"/info{0}.json".format(region),"w") as outfile:
        json.dump(info, outfile, indent=4, sort_keys=True)

exit()

# for prefire corrections

# save corrected MergedCSVperLS
print("store data - Muon prefiring")
data = data[['zDelBB_mc_prefMuon_nominal','zDelBE_mc_prefMuon_nominal','zDelEE_mc_prefMuon_nominal',
    'zDelBB_mc_prefECAL_nominal','zDelBE_mc_prefECAL_nominal','zDelEE_mc_prefECAL_nominal',
    u'recorded(/pb)', 'ls', 'time', 'run', 'measurement', 'fill', 'pileUp']]

data['zDelBB_mc'] = data['zDelBB_mc_prefMuon_nominal']
data['zDelBE_mc'] = data['zDelBE_mc_prefMuon_nominal']
data['zDelEE_mc'] = data['zDelEE_mc_prefMuon_nominal']
# data['zDelI_mc'] = data['zDelI_mc_prefMuon_nominal']

data_store = data[['ls', u'time', u'run', u'measurement', u'fill', u'pileUp',
       # u'delivered(/pb)',
       u'recorded(/pb)',
       # u'effBB_mc', u'RecBB', u'zRecBB', u'zDelBB',
       u'zDelBB_mc',
       # u'effBE_mc', u'RecBE', u'zRecBE', u'zDelBE',
       u'zDelBE_mc',
       # u'effEE_mc', u'RecEE', u'zRecEE', u'zDelEE',
       u'zDelEE_mc',
       # u'zDelI_mc',
       ]]

with open(outDir + '/Mergedcsvfile_perLS_corrMuonPrefire.csv', 'w') as file:
    data.to_csv(file, index=False)

print("store data - ECAL prefiring")

data['zDelBB_mc'] = data['zDelBB_mc_prefECAL_nominal']
data['zDelBE_mc'] = data['zDelBE_mc_prefECAL_nominal']
data['zDelEE_mc'] = data['zDelEE_mc_prefECAL_nominal']
# data['zDelI_mc'] = data['zDelI_mc_prefECAL_nominal']

data_store = data[['ls', u'time', u'run', u'measurement', u'fill', u'pileUp',
       # u'delivered(/pb)',
       u'recorded(/pb)',
       # u'effBB_mc', u'RecBB', u'zRecBB', u'zDelBB',
       u'zDelBB_mc',
       # u'effBE_mc', u'RecBE', u'zRecBE', u'zDelBE',
       u'zDelBE_mc',
       # u'effEE_mc', u'RecEE', u'zRecEE', u'zDelEE',
       u'zDelEE_mc',
       # u'zDelI_mc',       
       ]]

with open(outDir + '/Mergedcsvfile_perLS_corrECALPrefire.csv', 'w') as file:
    data.to_csv(file, index=False)

print("store data - done")
