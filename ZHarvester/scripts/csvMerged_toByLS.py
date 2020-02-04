import pandas as pd
import argparse
import pdb

parser = argparse.ArgumentParser()
parser.add_argument("csvFile", help="csvFile from ZCounting.py", type=str)
parser.add_argument("outputDir", help="outputDir", type=str)
args = parser.parse_args()


sigma_fid = 610.1401700042975
secPerLS = float(23.3)

ZBBAcc = 0.077904  # Fraction of Z events where both muons are in barrel region
ZBEAcc = 0.117200  # barrel and endcap
ZEEAcc = 0.105541  # both endcap

fullAcc = (ZBBAcc + ZBEAcc + ZEEAcc)

ZBBAcc = ZBBAcc / fullAcc
ZBEAcc = ZBEAcc / fullAcc
ZEEAcc = ZEEAcc / fullAcc

df = pd.read_csv(args.csvFile)

for key, source, factor in (('zDel_mc','ZCounts_MC', 1.),
                    ('zDelBB_mc','ZCountsBB_MC', ZBBAcc),
                    ('zDelEE_mc','ZCountsEE_MC', ZEEAcc),
                    ('zDel','ZCounts', 1.),
                    ('zDelBB','ZCountsBB', ZBBAcc),
                    ('zDelEE','ZCountsEE', ZEEAcc)):
    print(">>> write byLS_{0}.csv".format(source))

    df['delivered(hz/ub)'] = df[key] / (sigma_fid * secPerLS * factor) * 1000000

    with open(args.outputDir + '/byLS_{0}.csv'.format(source), 'wb') as file:
        file.write("#run:fill,ls,time,recorded(hz/ub),source\n")
        for run, fill, ls, ttime, ZinstLumiRec in df[['run','fill','ls','time','delivered(hz/ub)']].values:
            file.write("{0}:{1},{2}:{2},{3},{4},{5}\n".format(int(run), int(fill), int(ls), int(ttime), ZinstLumiRec, source))
