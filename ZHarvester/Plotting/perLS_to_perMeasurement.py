import os
import ROOT
import argparse
import pandas as pd
import uncertainties as unc
import pdb


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

from common import parsing
from common.logging import child_logger
log = child_logger(__name__)

parser = parsing.parser_plot()
parser.add_argument("-i", "--input", required=True, type=str, help="Nominator csv file with z rates per lumisection")
args = parser.parse_args()

outDir = args.output
if not os.path.isdir(outDir):
    os.mkdir(outDir)


########## Data Acquisition ##########
print("get Z luminosity")
data = pd.read_csv(args.input, sep=',',low_memory=False)

data['zDelBB_mc'] = data['zDelBB_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDelBE_mc'] = data['zDelBE_mc'].apply(lambda x: unc.ufloat_fromstr(x))
data['zDelEE_mc'] = data['zDelEE_mc'].apply(lambda x: unc.ufloat_fromstr(x))


data['zDel_mc'] = data['zDelBB_mc'] + data['zDelBE_mc'] + data['zDelEE_mc']


# data['zDelBB'] = data['zDelBB'].apply(lambda x: unc.ufloat_fromstr(x))
# data['zDelBE'] = data['zDelBE'].apply(lambda x: unc.ufloat_fromstr(x))
# data['zDelEE'] = data['zDelEE'].apply(lambda x: unc.ufloat_fromstr(x))
#
# data['zDel'] = data['zDelBB'] + data['zDelBE'] + data['zDelEE']

secPerLS = float(23.3)

print("make new perMeasurement csv")
results = []
for nFill, dFill in data.groupby("fill"):

    for nRun, dRun in dFill.groupby("run"):

        for nM, data_m in dRun.groupby("measurement"):

            # --- per measurement result
            # delLumi_m = data_m['delivered(/pb)'].sum()
            recLumi_m = data_m['recorded(/pb)'].sum()
            # deadtime_m = recLumi_m / delLumi_m
            timeWindow_m = len(data_m) * secPerLS

            res = {
                "fill": int(nFill),
                "run": int(nRun),
                "tdate_begin": min(data_m['time']),
                "tdate_end": max(data_m['time']),
                # "yieldBB": data_m['yieldBB'].sum(),
                # "yieldBE": data_m['yieldBE'].sum(),
                # "yieldEE": data_m['yieldEE'].sum(),
                # "zYieldBB": data_m['zYieldBB'].sum(),
                # "zYieldBE": data_m['zYieldBE'].sum(),
                # "zYieldEE": data_m['zYieldEE'].sum(),
                # "zDelBB": data_m['zDelBB'].sum(),
                # "zDelBE": data_m['zDelBE'].sum(),
                # "zDelEE": data_m['zDelEE'].sum(),
                "zDelBB_mc": data_m['zDelBB_mc'].sum(),
                "zDelBE_mc": data_m['zDelBE_mc'].sum(),
                "zDelEE_mc": data_m['zDelEE_mc'].sum(),
                # "xsecBB": data_m['zDelBB'].sum() / recLumi_m,
                # "xsecBE": data_m['zDelBE'].sum() / recLumi_m,
                # "xsecEE": data_m['zDelEE'].sum() / recLumi_m,
                # "xsecBB_mc": data_m['zDelBB_mc'].sum() / recLumi_m,
                # "xsecBE_mc": data_m['zDelBE_mc'].sum() / recLumi_m,
                # "xsecEE_mc": data_m['zDelEE_mc'].sum() / recLumi_m,
                # "lumiDel": delLumi_m,
                "lumiRec": recLumi_m,
                "timewindow": timeWindow_m,
                # "deadtime": deadtime_m,
                "pileUp": data_m['avgpu'].mean(),
                # "ZBBeff_mc": data_m['effBB_mc'].values.mean(),
                # "ZBEeff_mc": data_m['effBE_mc'].values.mean(),
                # "ZEEeff_mc": data_m['effEE_mc'].values.mean(),
            }

            results.append(res)

print("store data - done")

## Write per measurement csv file - one per run
results = pd.concat([pd.DataFrame([result]) for result in results])

outfilename = args.input.split("/")[-1].replace("perLS","perMeasurement")
with open(outDir + "/" + outfilename, 'w') as file:
    results.to_csv(file, index=False)

print("store data - done")
