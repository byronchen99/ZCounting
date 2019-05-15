import ROOT
import os,sys
from array import array
from operator import truediv
import pandas as pd
import numpy as np
import glob
import logging as log
import argparse
import pdb

parser = argparse.ArgumentParser()
parser.add_argument("-b","--beginRun",help="first run to analyze [%default]",default=299918)
parser.add_argument("-e","--endRun",help="analyze stops when comes to this run [%default]",default=1000000)
parser.add_argument("-p","--parametrizeType",help="define parametrization: 1 is for extrapolation, 2 is for piece-wise function",default=1)
parser.add_argument("-v","--verbose",help="increase logging level from INFO to DEBUG",default=False,action="store_true")
parser.add_argument("-c","--writeSummaryCSV",help="produce merged CSV with all runs",default=True)
parser.add_argument("-d","--dirDQM",help="Directory to the input root files from the DQM Offline module",default="/eos/home-d/dwalter/www/ZCounting/DQM-Offline-2018/")
parser.add_argument("-f","--ByLsCSV",help="ByLs csv input generated by testBril.sh",default="/eos/cms/store/group/comm_luminosity/ZCounting/brilcalcFile2018/*csv")
parser.add_argument("-g","--dirMC",help="Directory to root files for pilupe reweighting",default="/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/LookupTable/")
parser.add_argument("-t","--dirMCShape",help="Directory to root file for Z mass template",default="MCFiles/92X_norw_IsoMu27_noIso/")
parser.add_argument("-a","--dirCSV",help="where to write/store the CSV files",default="./")
parser.add_argument("-x","--dirEff",help="where to write/store efficiency Plots",default="./")

args = parser.parse_args()
if args.verbose:
    log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
else:
    log.basicConfig(format="%(levelname)s: %(message)s", level=log.INFO)

########## Input configuration ##########

#ByLS csv inputs generated by testBRIL.sh 
#inFile="/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/CloneJob/2017LumiByLS_hfet_trig_PU.csv"
inFileList=glob.glob(args.ByLsCSV)
inFileList.sort(key=os.path.getmtime)
inFile=inFileList[-1]
print "The brilcalc csv file: "+str(inFile)

eosDir= args.dirDQM 

#MC inputs: to build MC*Gaussian template for efficiency fitting

mcDir= args.dirMC # "/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/LookupTable/"
mcShapeSubDir= args.dirMCShape# "MCFiles/92X_norw_IsoMu27_noIso/"

########################################

########### Constant settings ##########
secPerLS=float(23.3)
currentYear=2018
maximumLS=2500
LSperMeasurement=100 #required number of lumi sections per measurement
staFitChi2Th=2.      #threshold on chi2 to trigger protection mechanism
staFitEffThHi=0.999  #threshold on eff. to trigger protection mechanism
staFitEffThLo=0.95   #threshold on eff. to trigger protection mechanism

ZfpRate = 0.01		#False positive rate of Z: background events counted as Z
ZBBRate = 0.077904 	#Fraction of Z events where both muons are in barrel region
ZBERate = 0.117200	# barrel and endcap
ZEERate = 0.105541	# both endcap

########################################

########## Pileup Correctins ###########

paraType=int(args.parametrizeType)

def fGen(a1,b1,a2,b2,a3,b3,fstep,para):
    #This function generates the pileup correction function, there are two models implemented
    #  para=1: linear function 
    #  para=2: stepwise linear function
    if para==1:
        def f(x):
            return a1 + b1*x
    elif para==2:
        def f(x):
            if x < fstep:
                return a2 + b2 * x
            else:
                return a3 + b3 * x
    else:
        print("ERROR: invalid parameterization type for pilup correction function")
    return f

#Define Z efficiency correction for pileup - parameters from MC
f_ZBBEffCorr = fGen(0.00305155,0.000519427, 0.0123732 ,0.000161345, -0.0878557 ,0.00218727 , 50, paraType)
f_ZBEEffCorr = fGen(0.00585766,0.000532229, 0.00875762,0.000493846, -0.0600895 ,0.00179    , 55, paraType)
f_ZEEEffCorr = fGen(0.0114659 ,0.00048351 , 0.0160629 ,0.000296308,  0.00132398,0.000743008, 40, paraType)

########################################

#log.info("Loading C marco...")	#I think we don't need this
ROOT.gROOT.LoadMacro(os.path.dirname(os.path.realpath(__file__))+"/calculateDataEfficiency.C") #load function getZyield(...) and calculateDataEfficiency(...)

#turn off graphical output on screen
ROOT.gROOT.SetBatch(True)

log.info("Loading input byls csv...")
lumiFile=open(str(inFile))
lumiLines=lumiFile.readlines()
data=pd.read_csv(inFile, sep=',',low_memory=False, skiprows=lambda x: lumiLines[x].startswith('#') and not lumiLines[x].startswith('#run'))
log.debug("%s",data.axes)
log.info("Loading input byls csv DONE...")
#formatting the csv

data[['run','fill']] = data['#run:fill'].str.split(':',expand=True).apply(pd.to_numeric)
data['ls'] = data['ls'].str.split(':',expand=True)[0].apply(pd.to_numeric)
data = data.drop(['#run:fill','hltpath','source'],axis=1)

if 'delivered(/ub)' in data.columns.tolist():      #convert to /pb
    data['delivered(/ub)'] = data['delivered(/ub)'].apply(lambda x:x / 1000000.)
    data['recorded(/ub)'] = data['recorded(/ub)'].apply(lambda x:x / 1000000.)
    data = data.rename(index=str, columns={'delivered(/ub)':'delivered(/pb)', 'recorded(/ub)':'recorded(/pb)' })

#if there are multiple entries of the same ls (for example from different triggers), only keep the one with the highest luminosity. 
data = data.sort_values(['fill','run','ls','delivered(/pb)','recorded(/pb)'])
data = data.drop_duplicates(['fill','run','ls'])

log.info("Looping over runs...")
for run in data.drop_duplicates('run')['run'].values:

    data_run = data.loc[data['run'] == run]

    fill = data_run.drop_duplicates('fill')['fill'].values[0]
    LSlist = data_run['ls'].values.tolist()

    if run<int(args.beginRun) or run>=int(args.endRun):
        continue
    
    #check if run was processed already
    processedRun = glob.glob(args.dirEff+'Run'+str(run))
    if len(processedRun)>0:
        print "Run "+str(run)+" was already processed, skipping and going to next run"
        continue

    #era split follows here:https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2017Analysis#DATA

    log.info("===Running Run %i",run)
    log.info("===Running Fill %i",fill)
        
    log.debug("===Setting up arrays for output csv...")
    fillarray=array('d')
    tdate_begin=[]
    tdate_end=[]
    Zrate=array('d')
    Zrate_EStatUp=array('d')
    Zrate_EStatDown=array('d')
    ZrateUncorrected=array('d')
    instDel=array('d')
    lumiDel=array('d')
    pileUp=array('d')
    ZyieldDel=array('d')

    ZyieldRec=array('d')
    lumiRec=array('d')
    windowarray=array('d')
    deadTime=array('d')
    beginLS=array('i')
    endLS=array('i')

    Zyield_fpr=array('d')    # Z fals positive rate: fraction of bkg events in reconstructed Zs: bkg/(sig+bkg)
    Zyield_chi2=array('d')

    # Efficiency related
    HLTeffB=array('d')
    HLTeffE=array('d')
    SITeffB=array('d')
    SITeffE=array('d')
    GloeffB=array('d')
    GloeffE=array('d')
    StaeffB=array('d')
    StaeffE=array('d')
    TrkeffB=array('d')
    TrkeffE=array('d')

    HLTeffB_chi2pass=array('d')
    HLTeffB_chi2fail=array('d')
    HLTeffE_chi2pass=array('d')
    HLTeffE_chi2fail=array('d')
    SITeffB_chi2pass=array('d')
    SITeffB_chi2fail=array('d')
    SITeffE_chi2pass=array('d')
    SITeffE_chi2fail=array('d')
    GloeffB_chi2pass=array('d')
    GloeffB_chi2fail=array('d')
    GloeffE_chi2pass=array('d')
    GloeffE_chi2fail=array('d')
    StaeffB_chi2pass=array('d')
    StaeffB_chi2fail=array('d')
    StaeffE_chi2pass=array('d')
    StaeffE_chi2fail=array('d')
    TrkeffB_chi2pass=array('d')
    TrkeffB_chi2fail=array('d')
    TrkeffE_chi2pass=array('d')
    TrkeffE_chi2fail=array('d')

    ZMCeff=array('d')
    ZMCeffBB=array('d')
    ZMCeffBE=array('d')
    ZMCeffEE=array('d')

    Zeff=array('d')
    ZBBeff=array('d')
    ZBEeff=array('d')
    ZEEeff=array('d')

    nMeasurements=0
    prevStaEffB=0.98
    prevStaEffE=0.98

    log.info("===Loading input DQMIO.root file...")
    eosFileList = glob.glob(eosDir+'/*/*'+str(run)+'*root')

    if not len(eosFileList)>0:
	print "The file does not yet exist for run: "+str(run)
	continue
    else:
	eosFile=eosFileList[0]

    print "The file exists: "+str(eosFile)+" for run  "+str(run)
    log.info("===Looping over measurements...")

    while len(LSlist) > 0: #begin next measurement "m"
        log.debug("Openning DQMIO.root file: %s", eosFile)

        f1 = ROOT.TFile(eosFile)
        h_yield_Z = f1.Get("DQMData/Run "+str(run)+"/ZCounting/Run summary/Histograms/h_mass_yield_Z").ProjectionX()
        # number of events in each ls (which may or may not have a Z candidate)
        h_n0 = f1.Get("DQMData/Run "+str(run)+"/ZCounting/Run summary/Histograms/h_npv").ProjectionX() 
        # produce goodLSlist with ls that are used for one measurement
        n0list = []
        goodLSlist = []
        while len(goodLSlist) < LSperMeasurement and len(LSlist) > 0:
            if len(LSlist) < 1:
                print("No more lumi sections in current run")
                break
            if LSlist[0] > h_yield_Z.GetNbinsX():
                log.error("======Lumi Section not stored in root file")
                break
            n0_ls = h_n0.GetBinContent(LSlist[0]+1)
            if n0_ls > 0:
                n0list.append(n0_ls)
                goodLSlist.append(LSlist[0])
            del LSlist[0]

        avgpu_m = sum(data_run.loc[data_run['ls'].isin(goodLSlist)]['avgpu'].values * np.array(n0list))/sum(n0list)
        
        recLumi_m = sum(data_run.loc[data_run['ls'].isin(goodLSlist)]['recorded(/pb)'].values)
        delLumi_m = sum(data_run.loc[data_run['ls'].isin(goodLSlist)]['delivered(/pb)'].values)
        deadtime_m = recLumi_m/delLumi_m
        timeWindow_m = len(goodLSlist) * secPerLS

        datestampLow_m = data_run.loc[data_run['ls'] == goodLSlist[0]]['time'].values[0].split(" ")
        datestampUp_m = data_run.loc[data_run['ls'] == goodLSlist[-1]]['time'].values[0].split(" ")

	dateLow_m=ROOT.TDatime(currentYear,int(datestampLow_m[0].split("/")[0]),int(datestampLow_m[0].split("/")[1]),int(datestampLow_m[1].split(":")[0]),int(datestampLow_m[1].split(":")[1]),int(datestampLow_m[1].split(":")[2]))
	dateUp_m=ROOT.TDatime(currentYear,int(datestampUp_m[0].split("/")[0]),int(datestampUp_m[0].split("/")[1]),int(datestampUp_m[1].split(":")[0]),int(datestampUp_m[1].split(":")[1]),int(datestampUp_m[1].split(":")[2]))
        
        log.debug("======beginTime: %s",dateLow_m.Convert())
        log.debug("======endTime: %s",dateUp_m.Convert())
        log.debug("======timeWindow: %f",timeWindow_m)

        Zyieldres_m = ROOT.getZyield(str(eosFile),args.dirEff,"h_mass_yield_Z",str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],1,5,recLumi_m);
        HLTeffresB_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"HLT",0,1,5,1,5,recLumi_m)
        HLTeffresE_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"HLT",1,1,5,1,5,recLumi_m)
        SITeffresB_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"SIT",0,1,1,1,1,recLumi_m)
        SITeffresE_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"SIT",1,1,1,1,1,recLumi_m)
        GloeffresB_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"Glo",0,2,2,2,2,recLumi_m,mcDir+mcShapeSubDir+"MuStaEff/MC/probes.root",mcDir)
        GloeffresE_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"Glo",1,2,2,2,2,recLumi_m,mcDir+mcShapeSubDir+"MuStaEff/MC/probes.root",mcDir)
        StaeffresB_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"Sta",0,2,2,2,2,recLumi_m,mcDir+mcShapeSubDir+"MuStaEff/MC/probes.root",mcDir)
        StaeffresE_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"Sta",1,2,2,2,2,recLumi_m,mcDir+mcShapeSubDir+"MuStaEff/MC/probes.root",mcDir)
        TrkeffresB_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"Trk",0,2,2,2,2,recLumi_m,mcDir+mcShapeSubDir+"MuStaEff/MC/probes.root",mcDir)
        TrkeffresE_m=ROOT.calculateDataEfficiency(str(eosFile),args.dirEff,str(run),nMeasurements,goodLSlist[0],goodLSlist[-1],avgpu_m,"Trk",1,2,2,2,2,recLumi_m,mcDir+mcShapeSubDir+"MuStaEff/MC/probes.root",mcDir)

        
        Zyield_m = Zyieldres_m[0]
        HLTeffB_m = HLTeffresB_m[0]
        HLTeffE_m = HLTeffresE_m[0]
        SITeffB_m = SITeffresB_m[0]
        SITeffE_m = SITeffresE_m[0]
        GloeffB_m = GloeffresB_m[0]
        GloeffE_m = GloeffresE_m[0]
        StaeffB_m = StaeffresB_m[0]
        StaeffE_m = StaeffresE_m[0]
        TrkeffB_m = TrkeffresB_m[0]
        TrkeffE_m = TrkeffresE_m[0]

        #if StaeffresB_m[3] > staFitChi2Th or StaeffresB_m[4] > staFitChi2Th or StaeffB_m >= staFitEffThHi or StaeffB_m <= staFitEffThLo:
        #    StaeffB_m = prevStaEffB
        #    log.warning("======Bad fit might happen, origin eff = %f, with chi2 = %f, %f",StaeffresB_m[0],StaeffresB_m[3],StaeffresB_m[4])
        #else:
        #    prevStaEffB = StaeffB_m

        #if StaeffresE_m[3] > staFitChi2Th or StaeffresE_m[4] > staFitChi2Th or StaeffE_m >= staFitEffThHi or StaeffE_m <= staFitEffThLo:
        #    StaeffE_m = prevStaEffE
        #    log.warning("======Bad fit might happen, origin eff = %f, with chi2 = %f, %f",StaeffresE_m[0],StaeffresE_m[3],StaeffresE_m[4])
        #else:
        #    prevStaEffE = StaeffE_m


        log.debug("======perMuonEff: %f, %f ,%f, %f, %f, %f",HLTeffB_m,HLTeffE_m,SITeffB_m,SITeffE_m,StaeffB_m,StaeffE_m)
        log.debug("======ZRawYield: %f",Zyield_m)

	#ZtoMuMu efficiency purely from data
        ZBBEff=(TrkeffB_m*TrkeffB_m * StaeffB_m*StaeffB_m * SITeffB_m*SITeffB_m * (1-(1-HLTeffB_m)*(1-HLTeffB_m)))
#        ZBBEff=(GloeffB_m*GloeffB_m * SITeffB_m*SITeffB_m * (1-(1-HLTeffB_m)*(1-HLTeffB_m)))
        ZBEEff=(TrkeffB_m*TrkeffE_m * StaeffB_m*StaeffE_m * SITeffB_m*SITeffE_m * (1-(1-HLTeffB_m)*(1-HLTeffE_m)))
#        ZBEEff=(GloeffB_m*GloeffE_m * SITeffB_m*SITeffE_m * (1-(1-HLTeffB_m)*(1-HLTeffE_m)))
        ZEEEff=(TrkeffE_m*TrkeffE_m * StaeffE_m*StaeffE_m * SITeffE_m*SITeffE_m * (1-(1-HLTeffE_m)*(1-HLTeffE_m)))
#        ZEEEff=(GloeffE_m*GloeffE_m * SITeffE_m*SITeffE_m * (1-(1-HLTeffE_m)*(1-HLTeffE_m)))

        #Statistic Uncertainties (low,high) error propagation 
        ZBBEff_EStat = [0.,0.]
        ZBEEff_EStat = [0.,0.]
        ZEEEff_EStat = [0.,0.]
        for i in (1,2):
            ZBBEff_EStat[i-1] = 2 * ZBBEff * np.sqrt( 
                            (TrkeffresB_m[i]/TrkeffB_m)**2 +
                            (StaeffresB_m[i]/StaeffB_m)**2 +
#                            (GloeffresB_m[i]/GloeffB_m)**2 +
                            (SITeffresB_m[i]/SITeffB_m)**2 + 
                            ((1-HLTeffB_m)/(1-(1-HLTeffB_m)**2)*HLTeffresB_m[i])**2 
                          )
            ZEEEff_EStat[i-1] = 2 * ZEEEff * np.sqrt(
                            (TrkeffresE_m[i]/TrkeffE_m)**2 +
                            (StaeffresE_m[i]/StaeffE_m)**2 +
#                            (GloeffresE_m[i]/GloeffE_m)**2 +
                            (SITeffresE_m[i]/SITeffE_m)**2 +
                            ((1-HLTeffE_m)/(1-(1-HLTeffE_m)**2)*HLTeffresE_m[i])**2 
                          )
            ZBEEff_EStat[i-1] = ZBEEff * np.sqrt( 
                            (TrkeffresB_m[i]/TrkeffB_m)**2 +                      
                            (TrkeffresE_m[i]/TrkeffE_m)**2 +
                            (StaeffresB_m[i]/StaeffB_m)**2 +
                            (StaeffresE_m[i]/StaeffE_m)**2 +
#                            (GloeffresB_m[i]/GloeffB_m)**2 +
#                            (GloeffresE_m[i]/GloeffE_m)**2 +
                            (SITeffresB_m[i]/SITeffB_m)**2 +
                            (SITeffresE_m[i]/SITeffE_m)**2 +
                            ((1-HLTeffE_m)/(1-(1-HLTeffB_m)*(1-HLTeffE_m))*HLTeffresB_m[i])**2 +
                            ((1-HLTeffB_m)/(1-(1-HLTeffB_m)*(1-HLTeffE_m))*HLTeffresE_m[i])**2 
                          )

	#ZtoMuMu efficiency correction as a parametrized function of pile-up
        ZBBEffCorr = f_ZBBEffCorr(avgpu_m)
	ZBEEffCorr = f_ZBEEffCorr(avgpu_m)
        ZEEEffCorr = f_ZEEEffCorr(avgpu_m)

	#ZtoMuMu efficiency after correction 
	ZMCEffBB = ZBBEff - ZBBEffCorr 
	ZMCEffBE = ZBEEff - ZBEEffCorr
	ZMCEffEE = ZEEEff - ZEEEffCorr

	#Multiply average frequency of each category with its efficiency
	ZMCEff = (ZMCEffBB * ZBBRate + ZMCEffBE * ZBERate + ZMCEffEE * ZEERate)/ (ZBBRate + ZBERate + ZEERate) 
        ZEff = (ZBBEff * ZBBRate + ZBEEff * ZBERate + ZEEEff * ZEERate)/ (ZBBRate + ZBERate + ZEERate)
        ZEff_EStat = [0.,0.]
        for i in (0,1):
            ZEff_EStat[i] = 1./(ZBBRate + ZBERate + ZEERate) * np.sqrt( 
                (ZBBRate*ZBBEff_EStat[i])**2 + (ZBERate*ZBEEff_EStat[i])**2 + (ZEERate*ZEEEff_EStat[i])**2 
                )

        #Or better take the actual frequency? -> bad statistic
        #ZMCEff = (ZMCEffBB*Zyield_BB_m + ZMCEffBE*(Zyield_m-Zyield_BB_m-Zyield_EE_m) + ZMCEffEE*Zyield_EE_m) / Zyield_m
        
        log.debug("======ZToMuMuEff: %f",ZMCEff)
        log.debug("======ZToMuMuEff: %f, %f ,%f, %f, %f, %f",ZMCEffBB, ZMCEffBE, ZMCEffEE, ZBBEff, ZBEEff, ZEEEff)

	#End products (about 1% fake rate)
#        ZMCXSec  = Zyield_m*(1-ZfpRate)/(ZMCEff*recLumi_m)
#        ZRateUncorrected  = Zyield_m*(1-ZfpRate)/(ZEff*timeWindow_m*deadtime_m)
#        ZRate  = Zyield_m*(1-ZfpRate)/(ZMCEff*timeWindow_m*deadtime_m)

        ZMCXSec  = Zyield_m/(ZMCEff*recLumi_m)
        ZRateUncorrected  = Zyield_m/(ZEff*timeWindow_m*deadtime_m)
        ZRate  = Zyield_m/(ZMCEff*timeWindow_m*deadtime_m)

        ZRate_EStat = [0.,0.]
        for i in (0,1):
            ZRate_EStat[i] = ZRate * ( Zyieldres_m[i]/Zyield_m + ZEff_EStat[i]/ZMCEff )
 
        log.debug("======ZMCXSec: %f",ZMCXSec)
        log.debug("======ZRate: %f",ZRate)
        
	#Variables to write in csv file
        fillarray.append(fill)

        tdate_begin.append(dateLow_m.Convert())
        tdate_end.append(dateUp_m.Convert())
        ZrateUncorrected.append(ZRateUncorrected)
        Zrate.append(ZRate)
        Zrate_EStatUp.append(ZRate_EStat[0])
        Zrate_EStatDown.append(ZRate_EStat[1])
        instDel.append(delLumi_m/timeWindow_m)
        lumiDel.append(delLumi_m)
	pileUp.append(avgpu_m)
        ZyieldDel.append(Zyield_m/(ZMCEff*deadtime_m))

	#Additional variables to write in efficiency csv file
        ZyieldRec.append(Zyield_m)
        lumiRec.append(recLumi_m)
        windowarray.append(timeWindow_m)
        deadTime.append(deadtime_m)
        beginLS.append(goodLSlist[0])
        endLS.append(goodLSlist[-1])

        Zyield_fpr.append(Zyieldres_m[4])
        Zyield_chi2.append(Zyieldres_m[3])

	#Efficiency related
    	HLTeffB.append(HLTeffB_m)
    	HLTeffE.append(HLTeffE_m)
        SITeffB.append(SITeffB_m)
        SITeffE.append(SITeffE_m)
        GloeffB.append(GloeffB_m)
        GloeffE.append(GloeffE_m)
        StaeffB.append(StaeffB_m)
        StaeffE.append(StaeffE_m)
        TrkeffB.append(TrkeffB_m)
        TrkeffE.append(TrkeffE_m)

        HLTeffB_chi2pass.append(HLTeffresB_m[3])
        HLTeffB_chi2fail.append(HLTeffresB_m[4])
        HLTeffE_chi2pass.append(HLTeffresE_m[3])
        HLTeffE_chi2fail.append(HLTeffresE_m[4])
        SITeffB_chi2pass.append(SITeffresB_m[3])
        SITeffB_chi2fail.append(SITeffresB_m[4])
        SITeffE_chi2pass.append(SITeffresE_m[3])
        SITeffE_chi2fail.append(SITeffresE_m[4])
        GloeffB_chi2pass.append(GloeffresB_m[3])
        GloeffB_chi2fail.append(GloeffresB_m[4])
        GloeffE_chi2pass.append(GloeffresE_m[3])
        GloeffE_chi2fail.append(GloeffresE_m[4])
        StaeffB_chi2pass.append(StaeffresB_m[3])
        StaeffB_chi2fail.append(StaeffresB_m[4])
        StaeffE_chi2pass.append(StaeffresE_m[3])
        StaeffE_chi2fail.append(StaeffresE_m[4])
        TrkeffB_chi2pass.append(TrkeffresB_m[3])
        TrkeffB_chi2fail.append(TrkeffresB_m[4])
        TrkeffE_chi2pass.append(TrkeffresE_m[3])
        TrkeffE_chi2fail.append(TrkeffresE_m[4])

        ZMCeff.append(ZMCEff)
        ZMCeffBB.append(ZMCEffBB)
        ZMCeffBE.append(ZMCEffBE)
        ZMCeffEE.append(ZMCEffEE)

        Zeff.append(ZEff)
        ZBBeff.append(ZBBEff)
        ZBEeff.append(ZBEEff)
        ZEEeff.append(ZEEEff)

        nMeasurements=nMeasurements+1


    ## Write Per Run CSV Files 
    print "Writing per Run CSV file"
    with open(args.dirCSV+'csvfile'+str(run)+'.csv','wb') as file:
        for c in range(0,nMeasurements):
		wline = (str(int(fillarray[c]))+","	+str(tdate_begin[c])+","	+str(tdate_end[c])+","
                        +str(Zrate[c])+","		+str(Zrate_EStatUp[c])+","	+str(Zrate_EStatDown[c])+","
                        +str(instDel[c])+","		+str(lumiDel[c])+","		+str(ZyieldDel[c])+","
                        +str(ZrateUncorrected[c]) )
                print(wline)
                file.write(wline + '\n')
    print "Writing per Run eff CSV file"
    with open(args.dirCSV+'effcsvfile'+str(run)+'.csv','wb') as file:
        for c in range(0,nMeasurements):
		wline = (str(int(fillarray[c]))+","	+str(tdate_begin[c])+"," 	+str(tdate_end[c])+","
                        +str(Zrate[c])+","         	+str(Zrate_EStatUp[c])+","	+str(Zrate_EStatDown[c])+","
                        +str(instDel[c])+","      	+str(lumiDel[c])+","   		+str(ZyieldDel[c])+","
                        +str(beginLS[c])+","       	+str(endLS[c])+","     		+str(lumiRec[c])+","  	+str(windowarray[c])+","
                        +str(Zyield_fpr[c])+","         +str(Zyield_chi2[c])+","
                        +str(HLTeffB[c])+","       	+str(HLTeffE[c])+","   
                        +str(SITeffB[c])+","       	+str(SITeffE[c])+","
                        +str(GloeffB[c])+","       	+str(GloeffE[c])+","
                        +str(StaeffB[c])+","            +str(StaeffE[c])+","
                        +str(TrkeffB[c])+","            +str(TrkeffE[c])+","
                        +str(HLTeffB_chi2pass[c])+","   +str(HLTeffB_chi2fail[c])+","
                        +str(HLTeffE_chi2pass[c])+","   +str(HLTeffE_chi2fail[c])+","
                        +str(SITeffB_chi2pass[c])+","	+str(SITeffB_chi2fail[c])+"," 
                        +str(SITeffE_chi2pass[c])+"," 	+str(SITeffE_chi2fail[c])+","
                        +str(GloeffB_chi2pass[c])+","	+str(GloeffB_chi2fail[c])+"," 
                        +str(GloeffE_chi2pass[c])+","  	+str(GloeffE_chi2fail[c])+","
                        +str(StaeffB_chi2pass[c])+","   +str(StaeffB_chi2fail[c])+","
                        +str(StaeffE_chi2pass[c])+","   +str(StaeffE_chi2fail[c])+","
                        +str(TrkeffB_chi2pass[c])+","   +str(TrkeffB_chi2fail[c])+","
                        +str(TrkeffE_chi2pass[c])+","   +str(TrkeffE_chi2fail[c])+","
                        +str(ZMCeff[c])+","        	+str(ZMCeffBB[c])+","  		+str(ZMCeffBE[c])+"," 	+str(ZMCeffEE[c])+","
                        +str(Zeff[c])+","               +str(ZBBeff[c])+","             +str(ZBEeff[c])+","    	+str(ZEEeff[c])+","
                        +str(pileUp[c]))
                file.write(wline + '\n')


## Write Big CSV File
print "Writing overall CSV file"
if args.writeSummaryCSV:
	rateFileList=sorted(glob.glob(args.dirCSV+'csvfile*.csv'))	
	with open(args.dirCSV+'Mergedcsvfile.csv','w') as file:
		file.write("fill,tdate_end,tdate_begin,ZRate,ZRate_EStatUp,ZRate_EStatDown,instDelLumi,delLumi,delZCount,ZRateUncorrected")
		file.write('\n')
		print "There are "+str(len(rateFileList))+" runs in the directory"
		for m in range(0,len(rateFileList)):						
			print "producing now csv file: "+rateFileList[m]
			try:
				iterF=open(rateFileList[m])
				lines=iterF.readlines()
				for line in lines:
					element=line.split(',')
					if element[3]=="nan" or element[3]=="0.0" or element[3]=="-0.0" or element[3]=="inf":
						continue
					file.write(line)
			except IOError as err:
    				print err.errno 
                                print err.strerror
				continue

        effFileList=sorted(glob.glob(args.dirCSV+'effcsvfile*.csv'))
	print "Starting to write efficiency files."	
        with open(args.dirCSV+'Mergedeffcsvfile.csv','wb') as fileTwo:
		fileTwo.write("fill,tdate_begin,tdate_end,ZRate,ZRate_EStatUp,ZRate_EStatDown,instDelLumi,delLumi,delZCount,beginLS,endLS,recLumi,windowarray,Zfpr,Zchi2,HLTeffB,HLTeffE,SITeffB,SITeffE,GloeffB,GloeffE,StaeffB,StaeffE,TrkeffB,TrkeffE,HLTeffB_chi2pass,HLTeffB_chi2fail,HLTeffE_chi2pass,HLTeffE_chi2fail,SITeffB_chi2pass,SITeffB_chi2fail,SITeffE_chi2pass,SITeffE_chi2fail,GloeffB_chi2pass,GloeffB_chi2fail,GloeffE_chi2pass,GloeffE_chi2fail,StaeffB_chi2pass,StaeffB_chi2fail,StaeffE_chi2pass,StaeffE_chi2fail,TrkeffB_chi2pass,TrkeffB_chi2fail,TrkeffE_chi2pass,TrkeffE_chi2fail,ZMCeff,ZMCeffBB,ZMCeffBE,ZMCeffEE,Zeff,ZBBeff,ZBEeff,ZEEeff,pileUp")
		fileTwo.write('\n')
		for m in range(0,len(effFileList)):

			print "producing now eff csv file: "+rateFileList[m]
			try:
				iterF=open(effFileList[m])
				lines=iterF.readlines()
				for line in lines:
					element=line.split(',')
					if element[3]=="nan" or element[3]=="0.0" or element[3]=="-0.0" or element[3]=="inf":
						continue
					fileTwo.write(line)
                        except IOError as err:
				print err.errno 
                                print err.strerror
                                continue
