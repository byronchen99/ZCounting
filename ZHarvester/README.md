# ZCounting

This code is to run over Z Counting histograms produced in DQM Offline and produce a CSV file with Z rates. 
The output file format has been agreed upon with ATLAS in 2016. 
It can be used to generate plots of the Z rate ratio of CMS/ATLAS. 
There are also python files to produce control plots of the Z reconstruction efficiencies and Z fiducial cross section.

A CMSSW environment needs to be set up to have the pandas library available.

## Z Rate Production
The ZCounting.py is the file to be executed to generate a CSV file with Z rates and additional information. 
Following parameters should be specified:

**--dirDQM DIRECTORY**\
Directory to the input root files from the DQM Offline module. 
Eg. the DQM Offline root files of the 2018 prompt-reco can be found [here](https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2018/SingleMuon/).

**--ByLsCSV BRILCALC.csv**\
A csv file with the information of run, fill, lumisection, time, delivered and recorded luminosity and avgpu is needed. The luminosity is used to produce control plots, the avgpu is needed for pileup corrections. A BRILCALC.csv file can be produced with 
```
./brilcalc\_makeByLs
```

**--dirCSV DIRECTORY**\
Specify a directory to store the output .csv file.

**--dirEff DIRECTORY**\
Specify a directory to store control plots of fits for the efficiency calculation.

**--dirMC DIRECTORY**\
The Sta efficiency is calculated by fitting a template function. The directory for the root files for generating this template have to be specified. For now, the default directory can be used.

**--dirMCShape DIRECTORY**\
Same as --dirMC

## Batch Processing
The script submit\_ZCounting can be used to run ZCounting.py on HTCondor.

## Plotting
The script cronMakePlots can be used to produce nice plots. The csv files from ZCounting.py output (and the ones from ATLAS) have to be specified. 
