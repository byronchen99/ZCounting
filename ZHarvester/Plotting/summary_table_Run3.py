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
# labels = "$\delta\NZ_{\mathrm{2022}}$"

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
info_binWidthUp   = load("/summary_nBins60/")
info_binWidthDown = load("/summary_nBins240/")


for key, iEra in info.items():
    # in case of empty info, continue
    if iEra['zDel'] == 0:
        continue
        
    keys.append(key)
    # nominal.append(iEra['prefECAL_nominal'])
    stat.append( 1./iEra['zDel'] )
    # stat.append((iEra['zDel_err']/iEra['zDel'] )**2
    #     + factor * (info['2017H']['zDel_err']/info['2017H']['zDel'])**2 )
    # bkgUp.append( (iEra['zRec_bkgUp']-iEra['zRec'])/iEra['zRec'] 
    #     - factor*(info['2017H']['zRec_bkgUp']-info['2017H']['zRec'])/info['2017H']['zRec'] )  # uncertainty on background subtraction
    # bkgDown.append( (iEra['zRec_bkgDown']-iEra['zRec'])/iEra['zRec'] 
    #     - factor*(info['2017H']['zRec_bkgDown']-info['2017H']['zRec'])/info['2017H']['zRec'] )

    cHLTUp.append( (iEra['zDel_cHLTUp']-iEra['zDel'])/iEra['zDel'] )  # uncertainty on HLT correlation
    cHLTDown.append( (iEra['zDel_cHLTDown']-iEra['zDel'])/iEra['zDel'] )
    cIOUp.append( (iEra['zDel_cIOUp']-iEra['zDel'])/iEra['zDel']  )  # uncertainty on HLT correlation
    cIODown.append( (iEra['zDel_cIODown']-iEra['zDel'])/iEra['zDel'] )

    # --- Muon prefire corrections       
    prefireSystUp.append(0.01)  # assume 1% uncertainty
    prefireSystDown.append(0.0)
    prefireECALUp.append(0.0)
    prefireECALDown.append(0.0)
    
    # --- from alternative sources
    # --- alternative signal
    # only take effect on efficiency
    # info_altSig1[key]['zDel'] = info_altSig1[key]['zDel'] / info_altSig1[key]['zRec'] * info[key]['zRec']
    # info_altSig2[key]['zDel'] = info_altSig2[key]['zDel'] / info_altSig2[key]['zRec'] * info[key]['zRec']
    if info_altSig1:
        altSig1.append((info_altSig1[key]['zDel']-info[key]['zDel'])/info[key]['zDel'])
    if info_altSig2:
        altSig2.append((info_altSig2[key]['zDel']-info[key]['zDel'])/info[key]['zDel'])
    if info_altBkg:                        
        altBkg.append((info_altBkg[key]['zDel']-info[key]['zDel'])/info[key]['zDel'])
    if info_lumiUp:
        lumiUp.append((info_lumiUp[key]['zDel']-info[key]['zDel'])/info[key]['zDel'])
    if info_lumiDown:
        lumiDown.append((info_lumiDown[key]['zDel']-info[key]['zDel'])/info[key]['zDel'])
    # if info_massUp:
        # massUp.append((info_massUp[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
            # -factor*)
    # if info_massDown:
        # massDown.append((info_massDown[key]['zDel']-info[key]['zDel'])/info[key]['zDel']
            # -factor*)
    if info_binWidthUp:
        binWidthUp.append((info_binWidthUp[key]['zDel']-info[key]['zDel'])/info[key]['zDel'])
    if info_binWidthDown:
        binWidthDown.append((info_binWidthDown[key]['zDel']-info[key]['zDel'])/info[key]['zDel'])

suffix = ""

outname = args.input.split("/")[-1].replace("info","result"+suffix)
outname = outname.replace(".json",".txt")
# ---- make latex table with variations:
with open(outDir+"/"+outname, "w") as outfile:

    columns = "l|"
    columns += "c"
    outfile.write(r"\begin{tabular}{"+columns+"}"+"\n")
    
    outfile.write(" & 2022" + r" \\"+"\n")

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
        ("Muon pref. up          ", prefireSystUp,   r" & $"),
        ("Muon pref. down        ", prefireSystDown, r" & $"),
        # ("ECAL pref. up          ", prefireECALUp,   r" & $"),
        # ("ECAL pref. down        ", prefireECALDown, r" & $"),
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
                entries.append(r"\NA $")
        outfile.write(name + connector + connector.join(entries) + r" \\"+"\n")
    outfile.write(r" \hline "+"\n")

    entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in stat]
    outfile.write(r"Statistical             & $\pm "+" & $\pm ".join(entries) + r" \\"+"\n")

    entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in sysUp]
    outfile.write(r"Systematic up           & $+"+" & $+".join(entries) + r" \\"+"\n")
    
    entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in sysDown]
    outfile.write(r"Systematic down         & $-"+" & $-".join(entries) + r" \\"+"\n")

    outfile.write(r" \hline "+"\n")
    outfile.write(r" \hline "+"\n")
    
    entries = [str(round(abs(np.sqrt(x + y))*100,2))+"$" for x, y in zip(stat, sysUp)]
    outfile.write(r"Total up                & $+"+" & $+".join(entries) + r" \\"+"\n")
    
    entries = [str(round(abs(np.sqrt(x + y))*100,2))+"$" for x, y in zip(stat, sysDown)]
    outfile.write(r"Total down              & $-"+" & $-".join(entries) + r" \\"+"\n")

    outfile.write(r"\end{tabular}"+"\n")
