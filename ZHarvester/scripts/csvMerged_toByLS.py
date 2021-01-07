import pandas as pd
import argparse
import pdb

parser = argparse.ArgumentParser()
parser.add_argument("csvFile", help="csvFile from ZCounting.py", type=str)
parser.add_argument("--xsec", help="csvFile from ZCounting.py where the cross section should be calculated from", type=str)

parser.add_argument("-o","--outputDir", help="outputDir", type=str)
args = parser.parse_args()

# --- get Z xsec
data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
xsecBB = sum(data_xsec['zDelBB_mc'])/sum(data_xsec['lumiRec'])
xsecBE = sum(data_xsec['zDelBE_mc'])/sum(data_xsec['lumiRec'])
xsecEE = sum(data_xsec['zDelEE_mc'])/sum(data_xsec['lumiRec'])

secPerLS = float(23.3)


df = pd.read_csv(args.csvFile)

for counts, source, xsec in (((df['zDelBB_mc']+df['zDelBE_mc']+df['zDelEE_mc']),'ZCounts', (xsecBB+xsecBE+xsecEE)),
                    ((df['zDelBB_mc']),'ZCountsBB', xsecBB),
                    ((df['zDelEE_mc']),'ZCountsEE', xsecEE)):
    print(">>> write byLS_{0}.csv".format(source))

    df['delivered(hz/ub)'] = counts / (xsec * secPerLS) * 1000000

    with open(args.outputDir + '/byLS_{0}.csv'.format(source), 'wb') as file:
        file.write("#dummy line\n")

        file.write("#run:fill,ls,time,recorded(hz/ub),source\n")
        for run, fill, ls, ttime, ZinstLumiRec in df[['run','fill','ls','time','delivered(hz/ub)']].values:
            if ZinstLumiRec == 0:
                print("No Z bosons in run {0}/ ls {1} ".format(run, ls))
                continue
            if int(run) < 305832:
                file.write("{0}:{1},{2}:{2},{3},{4},{5}\n".format(int(run), int(fill), int(ls), int(ttime), ZinstLumiRec, source))
            elif int(run) < 315257:
                file.write("{0}:{1},{2}:{2},{3},{4},{5}\n".format(int(run), int(fill), int(ls), int(ttime)-3600, ZinstLumiRec, source))
            else:   #for 2018 +1 year
                file.write("{0}:{1},{2}:{2},{3},{4},{5}\n".format(int(run), int(fill), int(ls), int(ttime)+31536000, ZinstLumiRec, source))

        file.write("#Summary:\n")
        file.write("#dummy line\n")
        file.write("#dummy line\n")
        file.write("#dummy line\n")
        file.write("#dummy line\n")
