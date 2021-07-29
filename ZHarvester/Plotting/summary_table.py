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
mc = []             # pileup dependent MC correction
prefireStatUp = []
prefireStatDown = []
prefireSystUp = []
prefireSystDown = []
prefireECALUp = []
prefireECALDown = []

for era, iEra in info.items():
    nominal.append(iEra['prefECAL_nominal'])
    stat.append(1./np.sqrt(iEra['yield']))
    mc.append((iEra['mc']-iEra['raw'])/iEra['mc'])
    prefireStatUp.append((iEra['prefMuon_StatUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])
    prefireStatDown.append((iEra['prefMuon_StatDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])
    prefireSystUp.append((iEra['prefMuon_SystUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])
    prefireSystDown.append((iEra['prefMuon_SystDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])
    prefireECALUp.append((iEra['prefECAL_Up']-iEra['prefECAL_nominal'])/iEra['prefECAL_nominal'])
    prefireECALDown.append((iEra['prefECAL_Down']-iEra['prefECAL_nominal'])/iEra['prefECAL_nominal'])
# ---- make latex table with variations:
with open(outDir+"/result.txt", "w") as outfile:

    n = len(info.keys())

    columns = "l|"
    columns += "".join(["c" for c in range(n)])
    outfile.write(r"\begin{tabular}{"+columns+"}"+"\n")

    outfile.write(" & " + " & ".join([x for x in info.keys()]) + r" \\"+"\n")

    outfile.write(r" \hline "+"\n")

    outfile.write("Nominal & "+" & ".join([str(abs(int(round(x,0)))) for x in nominal]) + r" \\" +"\n")

    outfile.write(r" \hline "+"\n")

    for name, values, connector in (
        ("Stat.                  ", stat, r" & $\pm "),
        ("MC corr.               ", mc, r" & $\pm "),
        ("muon pref. up (stat.)  ", prefireStatUp,   r" & $+"),
        ("muon pref. down (stat.)", prefireStatDown, r" & $"),
        ("muon pref. up (syst.)  ", prefireSystUp,   r" & $+"),
        ("muon pref. down (syst.)", prefireSystDown, r" & $"),
        ("ECAL pref. up          ", prefireECALUp,   r" & $+"),
        ("ECAL pref. down        ", prefireECALDown, r" & $")
    ):
        outfile.write(name + connector + connector.join([str(round(x*100,2))+"$" for x in values]) + r" \\"+"\n")

    outfile.write(r" \hline "+"\n")
    outfile.write(r"\end{tabular}"+"\n")
