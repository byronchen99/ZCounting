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
    "2017H": "$\delta\NZ_\lowPU$",
    "2018": "$\delta\NZ_{\mathrm{2018}}$"
},
{
    "2016preVFP": r"$\delta \frac{\NZ_\mathrm{2016 preVFP}}{\NZ_\lowPU}$ }",
    "2016postVFP": r"$\delta \frac{\NZ_\mathrm{2016 postVFP}}{\NZ_\lowPU}$ }",
    "2016H": r"$\delta \frac{\NZ_\mathrm{2016 H}}{\NZ_\lowPU}$ }",
    "2017": r"$\delta \frac{\NZ_\mathrm{2017}}{\NZ_\lowPU}$ }",
    "2018": r"$\delta \frac{\NZ_\mathrm{2018}}{\NZ_\lowPU}$ }"    
}]

# factor = 0 means we get total uncertainty on each period
# factor = 1 means we get uncertainty for the ratio to 2017H
for factor in (0, 1):

    keys = []
    nominal = []
    stat = []           # statistical uncertainty
    bkgUp = []          # background subtraction
    bkgDown = []
    cHLTUp = []             # pileup dependent di-muon HLT MC correction
    cHLTDown = []           # pileup dependent di-muon HLT MC correction
    cIOUp = []              # pileup dependent muon inner/outer track MC correction
    cIODown = []            # pileup dependent muon inner/outer track MC correction
    prefireStatUp = []
    prefireStatDown = []
    prefireSystUp = []
    prefireSystDown = []
    prefireECALUp = []
    prefireECALDown = []
    
    base = "/".join(args.input.split("/")[:-2])
    suffix = args.input.split("/")[-1]
    def load(name):
        if not os.path.isfile(base+name+suffix):
            return None
        with open(base+name+suffix, "r") as file_info:
            info = json.load(file_info)
        return info
        
    # ---  variations from alternative fits
    altSig1 = [] 
    altSig2 = []
    altBkg = []
    lumiUp = []
    lumiDown = []
    massUp = []
    massDown = []
    binWidthUp = []
    binWidthDown = []
    
    # load alternative infos
    info_altSig1      = load("/summary_altSigMCxCB/")
    info_altSig2      = load("/summary_altSigMC/")
    info_altBkg       = load("/summary_altBkg/")
    info_lumiUp       = load("/summary_lumi30/")
    info_lumiDown     = load("/summary_lumi15/")
    info_massUp       = load("/summary_massUp/")
    info_massDown     = load("/summary_massDown/")
    info_binWidthUp   = load("/summary_binWidth1/")
    info_binWidthDown = load("/summary_binWidth025/")

    # sortdict = {"2016preVFP":1, "2016postVFP":2, "2016H":3, "2017":4, "2017H":5, "2018":6}
    sortdict = {"2016preVFP":1, "2016postVFP":2, "2017":4, "2017H":5, "2018":6}
    
    items = filter(lambda (x, y): x in sortdict.keys(), info.items())

    # for alternative signal - only take the effect on the efficiency
    # info_altSig['2017H']['zDel'] = info_altSig['2017H']['zDel'] / info_altSig['2017H']['zRec'] * info['2017H']['zRec']
    
    for key, iEra in sorted(items, key=lambda item: sortdict[item[0]]):
        # in case of empty info, continue
        if iEra['zDel'] == 0:
            continue
        # if we compute the ratio with 2017H, we don't do it for 2017H :)
        if factor and key == "2017H":
            continue
            
        keys.append(key)
        # nominal.append(iEra['prefECAL_nominal'])
        stat.append( 1./iEra['zDel'] 
            + factor*1./info['2017H']['zDel'] )
        # stat.append((iEra['zDel_err']/iEra['zDel'] )**2
        #     + factor * (info['2017H']['zDel_err']/info['2017H']['zDel'])**2 )
        # bkgUp.append( (iEra['zRec_bkgUp']-iEra['zRec'])/iEra['zRec'] 
        #     - factor*(info['2017H']['zRec_bkgUp']-info['2017H']['zRec'])/info['2017H']['zRec'] )  # uncertainty on background subtraction
        # bkgDown.append( (iEra['zRec_bkgDown']-iEra['zRec'])/iEra['zRec'] 
        #     - factor*(info['2017H']['zRec_bkgDown']-info['2017H']['zRec'])/info['2017H']['zRec'] )

        cHLTUp.append( (iEra['zDel_cHLTUp']-iEra['zDel'])/iEra['zDel'] 
            - factor*(info['2017H']['zDel_cHLTUp']-info['2017H']['zDel'])/info['2017H']['zDel'] )  # uncertainty on HLT correlation
        cHLTDown.append( (iEra['zDel_cHLTDown']-iEra['zDel'])/iEra['zDel'] 
            - factor*(info['2017H']['zDel_cHLTDown']-info['2017H']['zDel'])/info['2017H']['zDel'] )
        cIOUp.append( (iEra['zDel_cIOUp']-iEra['zDel'])/iEra['zDel'] 
            - factor*(info['2017H']['zDel_cIOUp']-info['2017H']['zDel'])/info['2017H']['zDel'] )  # uncertainty on HLT correlation
        cIODown.append( (iEra['zDel_cIODown']-iEra['zDel'])/iEra['zDel'] 
            - factor*(info['2017H']['zDel_cIODown']-info['2017H']['zDel'])/info['2017H']['zDel'] )

        # --- Muon prefire corrections    
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
        
        # --- from alternative sources
        # --- alternative signal
        # only take effect on efficiency
        # info_altSig1[key]['zDel'] = info_altSig1[key]['zDel'] / info_altSig1[key]['zRec'] * info[key]['zRec']
        # info_altSig2[key]['zDel'] = info_altSig2[key]['zDel'] / info_altSig2[key]['zRec'] * info[key]['zRec']
        if info_altSig1:
            altSig1.append((info_altSig1[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                - factor*(info_altSig1['2017H']['zDel']-info['2017H']['zDel'])/info['2017H']['zDel'])
        if info_altSig2:
            altSig2.append((info_altSig2[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                - factor*(info_altSig2['2017H']['zDel']-info['2017H']['zDel'])/info['2017H']['zDel'])
        if info_altBkg:                        
            altBkg.append((info_altBkg[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                - factor*(info_altBkg['2017H']['zDel']-info['2017H']['zDel'])/info['2017H']['zDel'])
        if info_lumiUp:
            lumiUp.append((info_lumiUp[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                - factor*(info_lumiUp['2017H']['zDel']-info['2017H']['zDel'])/info['2017H']['zDel'])
        if info_lumiDown:
            lumiDown.append((info_lumiDown[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                - factor*(info_lumiDown['2017H']['zDel']-info['2017H']['zDel'])/info['2017H']['zDel'])
        # if info_massUp:
            # massUp.append((info_massUp[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                # -factor*)
        # if info_massDown:
            # massDown.append((info_massDown[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                # -factor*)
        if info_binWidthUp:
            binWidthUp.append((info_binWidthUp[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                - factor*(info_binWidthUp['2017H']['zDel']-info['2017H']['zDel'])/info['2017H']['zDel'])
        if info_binWidthDown:
            binWidthDown.append((info_binWidthDown[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
                - factor*(info_binWidthDown['2017H']['zDel']-info['2017H']['zDel'])/info['2017H']['zDel'])

    suffix = "_ratio" if factor==1 else ""

    outname = args.input.split("/")[-1].replace("info","result"+suffix)
    outname = outname.replace(".json",".txt")
    # ---- make latex table with variations:
    with open(outDir+"/"+outname, "w") as outfile:

        columns = "l|"
        columns += "".join(["c" for c in range(len(sortdict.keys()) - factor)])
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
        
        sysUp = [0 for x in stat]
        sysDown = [0 for x in stat]
        for name, values, connector in (
            ("Correlation \CHLT      ", cHLTUp,          r" & $\pm "),
            ("Correlation \CGlo      ", cIOUp,           r" & $\pm "),
            ("Background subtr.      ", bkgUp,   r" & $\pm "),
            ("Muon pref. up (stat.)  ", prefireStatUp,   r" & $"),
            ("Muon pref. down (stat.)", prefireStatDown, r" & $"),
            ("Muon pref. up (syst.)  ", prefireSystUp,   r" & $"),
            ("Muon pref. down (syst.)", prefireSystDown, r" & $"),
            ("ECAL pref. up          ", prefireECALUp,   r" & $"),
            ("ECAL pref. down        ", prefireECALDown, r" & $"),
            ("Alt. sig. model (MCxCB) ", altSig1,          r" & $"),
            ("Alt. sig. model (MC)    ", altSig2,          r" & $"),
            ("Alt. bkg. model        ", altBkg,          r" & $"),
            ("Lumi slice up          ", lumiUp,          r" & $"),
            ("Lumi slice down        ", lumiDown,        r" & $"),
            ("Mass range up          ", massUp,        r" & $"),
            ("Mass range down        ", massDown,      r" & $"),
            ("Bin width up           ", binWidthUp,      r" & $"),
            ("Bin width down         ", binWidthDown,    r" & $"),

        ):
            if len(values) == 0:
                continue
                
            entries = []
            for i, x in enumerate(values):
                if "\pm" in connector:
                    entries.append(str(round(abs(x)*100,2))+"$" )
                    sysUp[i] += x**2
                    sysDown[i] += x**2
                elif round(x*100,2) > 0:
                    entries.append("+"+str(round(x*100,2))+"$" )
                    sysUp[i] += x**2
                elif round(x*100,2) < 0:
                    entries.append(str(round(x*100,2))+"$" )
                    sysDown[i] += x**2
                else:
                    entries.append("\NA $")
            outfile.write(name + connector + connector.join(entries) + r" \\"+"\n")
        outfile.write(r" \hline "+"\n")

        entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in stat]
        outfile.write("Statistical             & $\pm "+" & $\pm ".join(entries) + r" \\"+"\n")

        entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in sysUp]
        outfile.write("Systematic up           & $+"+" & $+".join(entries) + r" \\"+"\n")
        
        entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in sysDown]
        outfile.write("Systematic down         & $-"+" & $-".join(entries) + r" \\"+"\n")

        outfile.write(r" \hline "+"\n")
        outfile.write(r" \hline "+"\n")
        
        entries = [str(round(abs(np.sqrt(x + y))*100,2))+"$" for x, y in zip(stat, sysUp)]
        outfile.write("Total up                & $+"+" & $+".join(entries) + r" \\"+"\n")
        
        entries = [str(round(abs(np.sqrt(x + y))*100,2))+"$" for x, y in zip(stat, sysDown)]
        outfile.write("Total down              & $-"+" & $-".join(entries) + r" \\"+"\n")

        outfile.write(r"\end{tabular}"+"\n")
