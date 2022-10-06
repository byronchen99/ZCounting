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
parser.add_argument("-x","--xsec",  default=None, type=str,
    help="csv file with z rates per lumisection where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Data Acquisition ##########
year = 2022
if year < 2022:
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


# data['zRec'] = data['zRecBB'] + data['zRecBE'] + data['zRecEE']
# data['zDelI'] = data['zDel'] #data['zDelBB'] + data['zDelBE'] + data['zDelEE']
# data['zDelUp'] = data['zDelBBUp'] + data['zDelBEUp'] + data['zDelEEUp']
# data['zDelDown'] = data['zDelBBDown'] + data['zDelBEDown'] + data['zDelEEDown']


data['effHLT']      = data['effHLT'].apply(lambda x: unc.ufloat_fromstr(x).n)
data['effSel']      = data['effSel'].apply(lambda x: unc.ufloat_fromstr(x).n)
data['effSta']      = data['effSta'].apply(lambda x: unc.ufloat_fromstr(x).n)
data['effGlo']      = data['effGlo'].apply(lambda x: unc.ufloat_fromstr(x).n)

# calculate correct statistical uncertainties - before HLT selection
data['NZ']       = data['zReco'].apply(lambda x: unc.ufloat_fromstr(x).n)
data['NbkgHLTFail'] = data['NbkgHLTFail'].apply(lambda x: unc.ufloat_fromstr(x).n)
data['NbkgHLTPass'] = data['NbkgHLTPass'].apply(lambda x: unc.ufloat_fromstr(x).n)

data['N1'] = 2*data['effHLT']*(1-data['cHLT']*data['effHLT'])*data['NZ'] + data['NbkgHLTFail']
data['N2'] = data['effHLT']**2 * data['cHLT'] * data['NZ'] + data['NbkgHLTPass']

data['N1'] = data['N1'].apply(lambda x: unc.ufloat(x, np.sqrt(x)) )
data['N2'] = data['N2'].apply(lambda x: unc.ufloat(x, np.sqrt(x)) )
data['N1bkg'] = data['NbkgHLTFail'].apply(lambda x: unc.ufloat(x, np.sqrt(x)) )
data['N2bkg'] = data['NbkgHLTPass'].apply(lambda x: unc.ufloat(x, np.sqrt(x)) )

data['N1sig'] = data['N1'] - data['N1bkg']
data['N2sig'] = data['N2'] - data['N2bkg']

data['effHLT'] = (2 * data['N2sig']) / (2 * data['N2sig'] + data['N1sig'])

# calculate correct statistical uncertainty for ID selection efficiency
for eff in ("Sel","Glo","Sta"):
    data['Nsig{0}'.format(eff)]     = data['Nsig{0}'.format(eff)].apply(lambda x: unc.ufloat_fromstr(x).n)
    data['Nbkg{0}Pass'.format(eff)] = data['Nbkg{0}Pass'.format(eff)].apply(lambda x: unc.ufloat_fromstr(x).n)
    data['Nbkg{0}Fail'.format(eff)] = data['Nbkg{0}Fail'.format(eff)].apply(lambda x: unc.ufloat_fromstr(x).n)

    data['N{0}Pass'.format(eff)] = data['Nsig{0}'.format(eff)]*data['eff{0}'.format(eff)] + data['Nbkg{0}Pass'.format(eff)]
    data['N{0}Fail'.format(eff)] = data['Nsig{0}'.format(eff)]*(1-data['eff{0}'.format(eff)]) + data['Nbkg{0}Fail'.format(eff)]

    data['N{0}Pass'.format(eff)] = data['N{0}Pass'.format(eff)].apply(lambda x: unc.ufloat(x, np.sqrt(x)) )
    data['N{0}Fail'.format(eff)] = data['N{0}Fail'.format(eff)].apply(lambda x: unc.ufloat(x, np.sqrt(x)) )
    data['Nbkg{0}Pass'.format(eff)] = data['Nbkg{0}Pass'.format(eff)].apply(lambda x: unc.ufloat(x, np.sqrt(x)) )
    data['Nbkg{0}Fail'.format(eff)] = data['Nbkg{0}Fail'.format(eff)].apply(lambda x: unc.ufloat(x, np.sqrt(x)) )

    data['Nsig{0}Pass'.format(eff)] = data['N{0}Pass'.format(eff)] - data['Nbkg{0}Pass'.format(eff)]
    data['Nsig{0}Fail'.format(eff)] = data['N{0}Fail'.format(eff)] - data['Nbkg{0}Fail'.format(eff)]

data['effSel'] = (2 * data['N2sig'] + data['N1sig']) / (2 * data['N2sig'] + data['N1sig'] + data['NsigSelFail'])
# data['effGlo'] = data["NsigSelPass"] / (data["NsigSelPass"]+data["NsigSelFail"])
# data['effGlo'] = data["NsigGloPass"] / (data["NsigGloPass"]+data["NsigGloFail"])
# data['effSta'] = data["NsigStaPass"] / (data["NsigStaPass"]+data["NsigStaFail"])

# ad-hoc uncertainties for background subtraction: vary background contribution by 50%
data['N1sig_Up'] = data['N1'] - 0.5*data['N1bkg']
data['N2sig_Up'] = data['N2'] - 0.5*data['N2bkg']

data['N1sig_Down'] = data['N1'] - 1.5*data['N1bkg']
data['N2sig_Down'] = data['N2'] - 1.5*data['N2bkg']

data['zRecI_bkgDown'] = (data['N2sig_Up'] + 0.5*data['N1sig_Up'])**2/data['N2sig_Up'] * data['cHLT']
data['zRecI_bkgDown'] = data['zRecI_bkgDown'].apply(lambda x: x.n)

data['zRecI_bkgUp'] = (data['N2sig_Down'] + 0.5*data['N1sig_Down'])**2/data['N2sig_Down'] * data['cHLT']
data['zRecI_bkgUp'] = data['zRecI_bkgUp'].apply(lambda x: x.n)

data['zRecI'] = (data['N2sig'] + 0.5*data['N1sig'])**2/data['N2sig'] * data['cHLT']
data['zDelI'] = (data['N2sig'] + 0.5*data['N1sig'])**2/data['N2sig'] * data['cHLT'] / (
    (2 * data['N2sig'] + data['N1sig']) / (2 * data['N2sig'] + data['N1sig'] + data['NsigSelFail'])*
    data['effSta']*data['effGlo']/data['cIO']
    )**2
data['zDelI_err'] = data['zDelI'].apply(lambda x: x.s)
data['zDelI'] = data['zDelI'].apply(lambda x: x.n)
data['zRecI'] = data['zRecI'].apply(lambda x: x.n)

# data['zRecI'] = data['zReco'].apply(lambda x: unc.ufloat_fromstr(x).n)
# data['zDelI'] = data['zDel'].apply(lambda x: unc.ufloat_fromstr(x).n)
# data['zDelI_err'] = data['zDel'].apply(lambda x: unc.ufloat_fromstr(x).s)

data['zDelI_cHLTUp'] = data['zDelI'] / data['cHLT'] * ((data['cHLT'] - 1)*1.5 + 1)
data['zDelI_cHLTDown'] = data['zDelI'] / data['cHLT'] * ((data['cHLT'] - 1)*0.5 + 1)

data['zDelI_cIOUp'] = data['zDelI'] / (data['cIO']**2) * ((data['cIO'] - 1)*2.0 + 1)**2
data['zDelI_cIODown'] = data['zDelI'] / (data['cIO']**2) * ((data['cIO'] - 1)*0.0 + 1)**2


regions = ("I",) #("","BB","BE","EE","I")
# --->>> prefire corrections
print("apply muon prefire corrections")
for var in prefire_variations_Muon:
    for region in regions:

        data['zDel{0}_prefMuon_{1}'.format(region,var)] = data['zDel{0}'.format(region)]

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
            data.loc[loc,'zDel{0}_prefMuon_{1}'.format(region,var)] *= factor

    # data['zDel_prefMuon_'+var] = data['zDelBB_prefMuon_'+var] \
    #     + data['zDelBE_prefMuon_'+var] \
    #     + data['zDelEE_prefMuon_'+var]

print("apply ECAL prefire corrections")

for var in prefire_variations_ECAL:
    for region in regions:

        data['zDel{0}_prefECAL_{1}'.format(region,var)] = data['zDel{0}_prefMuon_nominal'.format(region)]

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

            data.loc[loc, 'zDel{0}_prefECAL_{1}'.format(region,var)] *= func(data.loc[loc,'pileUp'], *params)

    # data['zDel_prefECAL_'+var] = data['zDelBB_prefECAL_'+var] \
    #     + data['zDelBE_prefECAL_'+var] \
    #     + data['zDelEE_prefECAL_'+var]

print("apply prefire corrections - done")
# <<<---

for region in regions:
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
        info[era]['zRec'] = sum(data.loc[loc,'zRec{0}'.format(region)])
        info[era]['zRec_bkgDown'] = sum(data.loc[loc,'zRec{0}_bkgDown'.format(region)])
        info[era]['zRec_bkgUp'] = sum(data.loc[loc,'zRec{0}_bkgUp'.format(region)])
        info[era]['zDel'] = sum(data.loc[loc,'zDel{0}'.format(region)])
        info[era]['zDel_err'] = np.sqrt(sum(data.loc[loc,'zDel{0}_err'.format(region)]**2))
        info[era]['zDel_cHLTUp']   = sum(data.loc[loc,'zDel{0}_cHLTUp'.format(region)])
        info[era]['zDel_cHLTDown'] = sum(data.loc[loc,'zDel{0}_cHLTDown'.format(region)])
        info[era]['zDel_cIOUp']    = sum(data.loc[loc,'zDel{0}_cIOUp'.format(region)])
        info[era]['zDel_cIODown']  = sum(data.loc[loc,'zDel{0}_cIODown'.format(region)])

        for var in prefire_variations_Muon:
            info[era]["prefMuon_"+var] = sum(data.loc[loc,'zDel{0}_prefMuon_{1}'.format(region, var)])
        for var in prefire_variations_ECAL:
            info[era]["prefECAL_"+var] = sum(data.loc[loc,'zDel{0}_prefECAL_{1}'.format(region, var)])

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
