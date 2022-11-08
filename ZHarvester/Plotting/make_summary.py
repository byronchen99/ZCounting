import os,sys
import argparse
import pandas as pd
import numpy as np
import uncertainties as unc
import pdb
import json

sys.path.append(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import pexp, pquad

parser = argparse.ArgumentParser()

parser.add_argument("-r","--rates", required=True, nargs='+', help="Nominator csv file with z rates per lumisection")
parser.add_argument("-x","--xsec",  default=None, type=str,
    help="csv file with z rates per lumisection where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Settings ##########

run2 = True
region = "I" #("","BB","BE","EE","I")

########## Data Acquisition ##########
if run2:
    # --- get prefire corrections
    prefire_variations_Muon = ('nominal', 'SystUp', 'SystDown', 'StatUp', 'StatDown')
    prefire_variations_ECAL = ("nominal", "Up", "Down",)
    with open('res/prefiring.json') as file_prefire:
        res_prefire = json.load(file_prefire)
else:
    prefire_variations_Muon = ()
    prefire_variations_ECAL = ()

# --- get Z xsec
if args.xsec:
    print("get Z cross section")
    data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)

    # --- z luminosity
    print("get Z luminosity")
    data = pd.concat([pd.read_csv(csv, sep=',',low_memory=False) for csv in args.rates] +[data_xsec,], ignore_index=True, sort=False)
else:
    # --- z luminosity
    print("get Z luminosity")
    data = pd.concat([pd.read_csv(csv, sep=',',low_memory=False) for csv in args.rates], ignore_index=True, sort=False)

data['recZCount_err'] = data['recZCount'].apply(lambda x: unc.ufloat_fromstr(x).s)
data['recZCount'] = data['recZCount'].apply(lambda x: unc.ufloat_fromstr(x).n)

data['recZCount_cHLTUp'] = data['recZCount'] / data['cHLT'] * ((data['cHLT'] - 1)*1.5 + 1)
data['recZCount_cHLTDown'] = data['recZCount'] / data['cHLT'] * ((data['cHLT'] - 1)*0.5 + 1)

data['recZCount_cIOUp'] = data['recZCount'] / (data['cIO']**2) * ((data['cIO'] - 1)*2.0 + 1)**2
data['recZCount_cIODown'] = data['recZCount'] / (data['cIO']**2) * ((data['cIO'] - 1)*0.0 + 1)**2

data['recZCount_cIDUp'] = data['recZCount'] / data['cID'] * ((data['cID'] - 1)*2.0 + 1)
data['recZCount_cIDDown'] = data['recZCount'] / data['cID'] * ((data['cID'] - 1)*0.0 + 1)

# number of selected Z bosons (before tracking efficiency corrections)
data['effTrk'] = data['effTrk'].apply(lambda x: unc.ufloat_fromstr(x).n)
data['selZCount'] = data['recZCount'] * data['effTrk']**2

# --->>> prefire corrections
print("apply muon prefire corrections")
for var in prefire_variations_Muon:

    data['recZCount_prefMuon_{0}'.format(var)] = data['recZCount']

    for lo, hi, era in (
        (272007, 278769, "2016preVFP"),
        (278769, 280919, "2016postVFP"),
        (280919, 284045, "2016H"),
        (297020, 306463, "2017"),
        (306828, 307083, "2017"),   # 2017H
        (315252, 325274, "2018"),
    ):
        loc = (data['run'] >= lo) & (data['run'] < hi)
        if len(loc) == 0:
            continue
            
        factor = float(res_prefire[era]["prefMuon"][region][var])
        data.loc[loc,'recZCount_prefMuon_{0}'.format(var)] *= factor


print("apply ECAL prefire corrections")

for var in prefire_variations_ECAL:

    data['recZCount_prefECAL_{0}'.format(var)] = data['recZCount']

    for lo, hi, era, func in (
        (272007, 278769, "2016preVFP", pexp),
        (278769, 284045, "2016postVFP", pexp),
        (297020, 306463, "2017", pexp),
        (306828, 307083, "2017H", pquad),   # 2017H
    ):
        loc = (data['run'] >= lo) & (data['run'] < hi)
    
        if len(loc) == 0:
            continue

        params = res_prefire[era]["prefECAL"][region][var]

        data.loc[loc, 'recZCount_prefECAL_{0}'.format(var)] *= func(data.loc[loc,'pileUp'], *params)

print("apply prefire corrections - done")
# <<<---

info = {}

for lo, hi, era in (
    (272007, 278769, "2016preVFP"),
    (278769, 284045, "2016postVFP"),
    (278769, 280919, "2016G"),
    (280919, 284045, "2016H"),
    (297020, 306463, "2017"),
    (306828, 307083, "2017H"),   # 2017H
    (315252, 325274, "2018"),
    (356100, 356616, "2022"),
):
    loc = (data['run'] >= lo) & (data['run'] < hi)

    if len(loc) == 0:
        continue

    info[era] = {}        
    info[era]['selZCount'] = sum(data.loc[loc,'selZCount'])
    info[era]['recZCount'] = sum(data.loc[loc,'recZCount'])
    info[era]['recZCount_err'] = np.sqrt(sum(data.loc[loc,'recZCount_err']**2))
    info[era]['recZCount_cHLTUp']   = sum(data.loc[loc,'recZCount_cHLTUp'])
    info[era]['recZCount_cHLTDown'] = sum(data.loc[loc,'recZCount_cHLTDown'])
    info[era]['recZCount_cIOUp']    = sum(data.loc[loc,'recZCount_cIOUp'])
    info[era]['recZCount_cIODown']  = sum(data.loc[loc,'recZCount_cIODown'])
    info[era]['recZCount_cIDUp']    = sum(data.loc[loc,'recZCount_cIDUp'])
    info[era]['recZCount_cIDDown']  = sum(data.loc[loc,'recZCount_cIDDown'])


    for var in prefire_variations_Muon:
        info[era]["prefMuon_"+var] = sum(data.loc[loc,'recZCount_prefMuon_{0}'.format(var)])
    for var in prefire_variations_ECAL:
        info[era]["prefECAL_"+var] = sum(data.loc[loc,'recZCount_prefECAL_{0}'.format(var)])

with open(outDir+"/info{0}.json".format(region),"w") as outfile:
    json.dump(info, outfile, indent=4, sort_keys=True)

exit()

# for prefire corrections

# save corrected MergedCSVperLS
print("store data - Muon prefiring")
data = data[['zDelBB_prefMuon_nominal','zDelBE_prefMuon_nominal','zDelEE_prefMuon_nominal',
    'zDelBB_prefECAL_nominal','zDelBE_prefECAL_nominal','zDelEE_prefECAL_nominal',
    u'recorded(/pb)', 'ls', 'time', 'run', 'measurement', 'fill', 'pileUp']]

data['zDelBB'] = data['zDelBB_prefMuon_nominal']
data['zDelBE'] = data['zDelBE_prefMuon_nominal']
data['zDelEE'] = data['zDelEE_prefMuon_nominal']
# data['zDelI'] = data['zDelI_prefMuon_nominal']

data_store = data[['ls', u'time', u'run', u'measurement', u'fill', u'pileUp',
       # u'delivered(/pb)',
       u'recorded(/pb)',
       # u'effBB', u'RecBB', u'zRecBB', u'zDelBB',
       u'zDelBB',
       # u'effBE', u'RecBE', u'zRecBE', u'zDelBE',
       u'zDelBE',
       # u'effEE', u'RecEE', u'zRecEE', u'zDelEE',
       u'zDelEE',
       # u'zDelI',
       ]]

with open(outDir + '/Mergedcsvfile_perLS_corrMuonPrefire.csv', 'w') as file:
    data.to_csv(file, index=False)

print("store data - ECAL prefiring")

data['zDelBB'] = data['zDelBB_prefECAL_nominal']
data['zDelBE'] = data['zDelBE_prefECAL_nominal']
data['zDelEE'] = data['zDelEE_prefECAL_nominal']
# data['zDelI'] = data['zDelI_prefECAL_nominal']

data_store = data[['ls', u'time', u'run', u'measurement', u'fill', u'pileUp',
       # u'delivered(/pb)',
       u'recorded(/pb)',
       # u'effBB', u'RecBB', u'zRecBB', u'zDelBB',
       u'zDelBB',
       # u'effBE', u'RecBE', u'zRecBE', u'zDelBE',
       u'zDelBE',
       # u'effEE', u'RecEE', u'zRecEE', u'zDelEE',
       u'zDelEE',
       # u'zDelI',       
       ]]

with open(outDir + '/Mergedcsvfile_perLS_corrECALPrefire.csv', 'w') as file:
    data.to_csv(file, index=False)

print("store data - done")
