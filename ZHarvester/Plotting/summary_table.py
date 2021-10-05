import json
import numpy as np
import argparse
import os
import pdb

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", required=True, type=str, help="input json from make_summary.py")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

with open(args.input, "r") as file_info:
    info = json.load(file_info)

nominal = []
stat = []           # statistical uncertainty
mcUp = []             # pileup dependent MC correction
mcDown = []           # pileup dependent MC correction
prefireStatUp = []
prefireStatDown = []
prefireSystUp = []
prefireSystDown = []
prefireECALUp = []
prefireECALDown = []

sortdict = {"2016preVFP":1, "2016postVFP":2, "2016H":3, "2017":4, "2017H":5, "2018":6}

for key, iEra in sorted(info.items(), key=lambda item: sortdict[item[0]]):
    nominal.append(iEra['prefECAL_nominal'])
    stat.append(1./np.sqrt(iEra['zRec_mc']))
    mcUp.append((iEra['mcUp']-iEra['mc'])/iEra['mc'])
    mcDown.append((iEra['mcDown']-iEra['mc'])/iEra['mc'])
    prefireStatUp.append((iEra['prefMuon_StatUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])
    prefireStatDown.append((iEra['prefMuon_StatDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])
    prefireSystUp.append((iEra['prefMuon_SystUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])
    prefireSystDown.append((iEra['prefMuon_SystDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])
    prefireECALUp.append((iEra['prefECAL_Up']-iEra['prefECAL_nominal'])/iEra['prefECAL_nominal'])
    prefireECALDown.append((iEra['prefECAL_Down']-iEra['prefECAL_nominal'])/iEra['prefECAL_nominal'])

outname = args.input.split("/")[-1].replace("info","result")
outname = outname.replace(".json",".txt")
# ---- make latex table with variations:
with open(outDir+"/"+outname, "w") as outfile:

    n = len(info.keys())

    columns = "l|"
    columns += "".join(["c" for c in range(n)])
    outfile.write(r"\begin{tabular}{"+columns+"}"+"\n")

    outfile.write(" & " + " & ".join([x for x, y in sorted(info.items(), key=lambda item: sortdict[item[0]])]) + r" \\"+"\n")

    outfile.write(r" \hline "+"\n")

    # outfile.write("Nominal & "+" & ".join([str(abs(int(round(x,0)))) for x in nominal]) + r" \\" +"\n")
    #
    # outfile.write(r" \hline "+"\n")

    for name, values, connector in (
        ("Statatistical          ", stat, r" & $\pm "),
        ("Correlation up         ", mcUp, r" & $+ "),
        ("Correlation down       ", mcDown, r" & $ "),
        ("Muon pref. up (stat.)  ", prefireStatUp,   r" & $+"),
        ("Muon pref. down (stat.)", prefireStatDown, r" & $"),
        ("Muon pref. up (syst.)  ", prefireSystUp,   r" & $+"),
        ("Muon pref. down (syst.)", prefireSystDown, r" & $"),
        ("ECAL pref. up          ", prefireECALUp,   r" & $+"),
        ("ECAL pref. down        ", prefireECALDown, r" & $")
    ):
        outfile.write(name + connector + connector.join([str(round(x*100,2))+"$" for x in values]) + r" \\"+"\n")

    outfile.write(r" \hline "+"\n")
    outfile.write(r"\end{tabular}"+"\n")
