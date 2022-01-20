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

# labels for table header
labels = [{
    "2016preVFP": "$\delta\NZ_{\mathrm{2016 preVFP}}$",
    "2016postVFP": "$\delta\NZ_{\mathrm{2016 postVFP}}$",
    "2016H": "$\delta\NZ_{\mathrm{2016 H}}$",
    "2017": "$\delta\NZ_{\mathrm{2017}}$",
    "2017H": "$\delta\NZ_{\mathrm{2017H}}$",
    "2018": "$\delta\NZ_{\mathrm{2018}}$"
},
{
    "2016preVFP": r"$\delta \frac{\NZ_\mathrm{2016 preVFP}}{\NZ_{\mathrm{2017H}}}$ }",
    "2016postVFP": r"$\delta \frac{\NZ_\mathrm{2016 postVFP}}{\NZ_{\mathrm{2017H}}}$ }",
    "2016H": r"$\delta \frac{\NZ_\mathrm{2016 H}}{\NZ_{\mathrm{2017H}}}$ }",
    "2017": r"$\delta \frac{\NZ_\mathrm{2017}}{\NZ_{\mathrm{2017H}}}$ }",
    "2018": r"$\delta \frac{\NZ_\mathrm{2018}}{\NZ_{\mathrm{2017H}}}$ }"    
}]

# factor = 0 means we get total uncertainty on each period
# factor = 1 means we get uncertainty for the ratio to 2017H
for factor in (0, 1):

    keys = []
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

    # variations from alternative fits
    altSig = []
    altBkg = []
    lumiUp = []
    lumiDown = []
    massUp = []
    massDown = []
    binWidthUp = []
    binWidthDown = []
    base = "/nfs/dust/cms/user/dwalter/data/Lumi/V13_05/"
    with open(base+"/summary_v2_altSig/"+args.input.split("/")[-1], "r") as file_info:
        info_altSig = json.load(file_info)
    with open(base+"/summary_v2_altBkg/"+args.input.split("/")[-1], "r") as file_info:
        info_altBkg = json.load(file_info)
    with open(base+"/summary_v2_lumiUp/"+args.input.split("/")[-1], "r") as file_info:
        info_lumiUp = json.load(file_info)
    with open(base+"/summary_v2_lumiDown/"+args.input.split("/")[-1], "r") as file_info:
        info_lumiDown = json.load(file_info)
    # with open(base+"/summary_massUp_v2/"+args.input.split("/")[-1], "r") as file_info:
    #     info_massUp = json.load(file_info)
    # with open(base+"/summary_massDown_v2/"+args.input.split("/")[-1], "r") as file_info:
    #     info_massDown = json.load(file_info)
    with open(base+"/summary_v2_binWidthUp/"+args.input.split("/")[-1], "r") as file_info:
        info_binWidthUp = json.load(file_info)
    with open(base+"/summary_v2_binWidthDown/"+args.input.split("/")[-1], "r") as file_info:
        info_binWidthDown = json.load(file_info)

    sortdict = {"2016preVFP":1, "2016postVFP":2, "2016H":3, "2017":4, "2017H":5, "2018":6}

    for key, iEra in sorted(info.items(), key=lambda item: sortdict[item[0]]):
        if factor and key == "2017H":
            continue
            
        keys.append(key)
        # nominal.append(iEra['prefECAL_nominal'])
        stat.append(np.sqrt(1./iEra['zRec_mc'] 
            + factor*1./info['2017H']['zRec_mc']) )
        mcUp.append( 0.5*(iEra['mcUp']-iEra['mc'])/iEra['mc'] 
            - factor*0.5*(info['2017H']['mcUp']-info['2017H']['mc'])/info['2017H']['mc'] )  # 50% uncertainty on HLT correlation
        # mcDown.append((iEra['mcDown']-iEra['mc'])/iEra['mc'] 
        #     - factor*(info['2017H']['mcDown']-info['2017H']['mc'])/info['2017H']['mc'])

        # Muon prefire corrections    
        # correlated between 2017 and 2017H because it's the same statistics used in the measurement of the prefire correction
        if key == "2017": 
            prefireStatUp.append((iEra['prefMuon_StatUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']
                - factor*((iEra['prefMuon_StatUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']) )
            prefireStatDown.append((iEra['prefMuon_StatDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']
                - factor*((iEra['prefMuon_StatDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']) )
        else:
            prefireStatUp.append(np.sqrt( ((iEra['prefMuon_StatUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])**2
                + factor*((iEra['prefMuon_StatUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])**2 ))
            prefireStatDown.append(-np.sqrt( ((iEra['prefMuon_StatDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])**2
                + factor*((iEra['prefMuon_StatDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])**2 ))
        
        prefireSystUp.append((iEra['prefMuon_SystUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']
            - factor*(info['2017H']['prefMuon_SystUp']-info['2017H']['prefMuon_nominal'])/info['2017H']['prefMuon_nominal']
        )
        prefireSystDown.append((iEra['prefMuon_SystDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']
            - factor*(info['2017H']['prefMuon_SystDown']-info['2017H']['prefMuon_nominal'])/info['2017H']['prefMuon_nominal']
        )
        prefireECALUp.append((iEra['prefECAL_Up']-iEra['prefECAL_nominal'])/iEra['prefECAL_nominal']
            - factor*(info['2017H']['prefECAL_Up']-info['2017H']['prefECAL_nominal'])/info['2017H']['prefECAL_nominal']
        )
        prefireECALDown.append((iEra['prefECAL_Down']-iEra['prefECAL_nominal'])/iEra['prefECAL_nominal']
            - factor*(info['2017H']['prefECAL_Down']-info['2017H']['prefECAL_nominal'])/info['2017H']['prefECAL_nominal']
        )
        altSig.append((info_altSig[key]['mc']-iEra['mc'])/iEra['mc']
            - factor*(info_altSig['2017H']['mc']-info['2017H']['mc'])/info['2017H']['mc'])
        altBkg.append((info_altBkg[key]['mc']-iEra['mc'])/iEra['mc']
            - factor*(info_altBkg['2017H']['mc']-info['2017H']['mc'])/info['2017H']['mc'])
        lumiUp.append((info_lumiUp[key]['mc']-iEra['mc'])/iEra['mc']
            - factor*(info_lumiUp['2017H']['mc']-info['2017H']['mc'])/info['2017H']['mc'])
        lumiDown.append((info_lumiDown[key]['mc']-iEra['mc'])/iEra['mc']
            - factor*(info_lumiDown['2017H']['mc']-info['2017H']['mc'])/info['2017H']['mc'])
        # massUp.append((info_massUp[key]['mc']-iEra['mc'])/iEra['mc']
            # -factor*)
        # massDown.append((info_massDown[key]['mc']-iEra['mc'])/iEra['mc']
            # -factor*)
        binWidthUp.append((info_binWidthUp[key]['mc']-iEra['mc'])/iEra['mc']
            - factor*(info_binWidthUp['2017H']['mc']-info['2017H']['mc'])/info['2017H']['mc'])
        binWidthDown.append((info_binWidthDown[key]['mc']-iEra['mc'])/iEra['mc']
            - factor*(info_binWidthDown['2017H']['mc']-info['2017H']['mc'])/info['2017H']['mc'])

    suffix = "_ratio" if factor==1 else ""

    outname = args.input.split("/")[-1].replace("info","result"+suffix)
    outname = outname.replace(".json",".txt")
    # ---- make latex table with variations:
    with open(outDir+"/"+outname, "w") as outfile:


        columns = "l|"
        columns += "".join(["c" for c in info.keys()])
        outfile.write(r"\begin{tabular}{"+columns+"}"+"\n")
        
        if factor:
            outfile.write(" \multirow{3}{*}{ } & \multirow{3}{*}{ " + " & \multirow{3}{*}{ ".join([labels[factor][x] for x in keys]) + r" \\"+"\n")
            outfile.write(" & " + " & ".join([" " for x in keys]) + r" \\"+"\n")
            outfile.write(" & " + " & ".join([" " for x in keys]) + r" \\"+"\n")
        else:
            outfile.write(" & " + " & ".join([labels[factor][x] for x in keys]) + r" \\"+"\n")

        outfile.write(r" \hline "+"\n")

        # outfile.write("Nominal & "+" & ".join([str(abs(int(round(x,0)))) for x in nominal]) + r" \\" +"\n")
        #
        # outfile.write(r" \hline "+"\n")
        
        totalsUp = [0 for x in stat]
        totalsDown = [0 for x in stat]
        for name, values, connector in (
            ("Statistical            ", stat, r" & $\pm "),
            ("HLT correlation        ", mcUp, r" & $\pm "),
            # ("HLT correlation       ", mcDown, r" & $ "),
            ("Alt. bkg. model        ", altBkg,   r" & $"),
            ("Alt. sig. model        ", altSig,   r" & $"),
            ("Lumi slice up          ", lumiUp,   r" & $"),
            ("Lumi slice down        ", lumiDown,   r" & $"),
            # ("Mass range up          ", massUp,   r" & $"),
            # ("Mass range down        ", massDown,   r" & $"),
            ("Bin width up           ", binWidthUp,   r" & $"),
            ("Bin width down         ", binWidthDown,   r" & $"),
            ("Muon pref. up (stat.)  ", prefireStatUp,   r" & $"),
            ("Muon pref. down (stat.)", prefireStatDown, r" & $"),
            ("Muon pref. up (syst.)  ", prefireSystUp,   r" & $"),
            ("Muon pref. down (syst.)", prefireSystDown, r" & $"),
            ("ECAL pref. up          ", prefireECALUp,   r" & $"),
            ("ECAL pref. down        ", prefireECALDown, r" & $")
        ):
            entries = []
            for i, x in enumerate(values):
                if "\pm" in connector:
                    entries.append(str(round(abs(x)*100,2))+"$" )
                    totalsUp[i] += x**2
                    totalsDown[i] += x**2
                elif round(x*100,2) > 0:
                    entries.append("+"+str(round(x*100,2))+"$" )
                    totalsUp[i] += x**2
                elif round(x*100,2) < 0:
                    entries.append(str(round(x*100,2))+"$" )
                    totalsDown[i] += x**2
                else:
                    entries.append("\NA $")
            outfile.write(name + connector + connector.join(entries) + r" \\"+"\n")
        outfile.write(r" \hline "+"\n")
        
        entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in totalsUp]
        outfile.write("TotalUp                 & $+"+" & $+".join(entries) + r" \\"+"\n")
        
        entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in totalsDown]
        outfile.write("TotalDown               & $-"+" & $-".join(entries) + r" \\"+"\n")

        outfile.write(r"\end{tabular}"+"\n")
