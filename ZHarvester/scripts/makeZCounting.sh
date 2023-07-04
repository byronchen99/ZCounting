#!/bin/bash

##############
# Set Variables
#############

pass=''                                                                                   # grid passwd
cmsswRel='/afs/cern.ch/user/'${USER:0:1}'/'$USER'/work/private/CMSSW_12_4_12/src/'             # path to ZCounting dir.
PDs=('/Muon0/' '/Muon1/')                                                                  # options for PD
COREDIR='https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2023/'
certif=/tmp/x509up_u81013                                                                         # you get one after the voms (grid proxy)
targetDir='/eos/cms/store/group/comm_luminosity/ZCounting/2023/DQMFiles/'
#output_ZCounting=''
output_ZCounting='/afs/cern.ch/user/'${USER:0:1}'/'$USER'/work/private/2023_Data_ZCounting'     # update to you directory 
normtagDate=$(date +'%Y%m%d')
normtagFileName="normtag_BRIL_${normtagDate}.json"
normtagFilePath="/eos/cms/store/group/comm_luminosity/ZCounting/2023/Normtags/$normtagFileName"

beginRun=366442
endRun=368823
etaMin=0.0
etaMax=2.4
brilcalcOutputFileName="byLS_Collisions23_${beginRun}_${endRun}_Muon_${normtagDate}.csv"
brilcalcOutputFilePath="/eos/cms/store/group/comm_luminosity/ZCounting/2023/brilcalcByLS/${brilcalcOutputFileName}"

##################
# Get List Of Files
#################

cd $cmsswRel
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
echo $pass | voms-proxy-init -voms cms -rfc    # get the grid proxy
cd /afs/cern.ch/user/"${USER:0:1}"/$USER/work/private/CMSSW_12_4_12/src/ZCounting/  # enter the ZCounting dir.
rm -rf filelist_*.log


# Create the log files with the files from the PDDIR (uses grid certs. and grid passwd)
for index in "${!PDs[@]}"
do
    PD=${PDs[$index]}
    PDDIR=$COREDIR$PD
    filelist="filelist_$index.log"
    for rundir in `curl -k --cert $certif -X GET $PDDIR | awk 's=index($0,"/000") { print substr($0,s+1,10)}'`
    do
        curl -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem --pass $pass -X GET $PDDIR$rundir | awk  -v COREDIR="$COREDIR" -F '<tr><td>' '{ print COREDIR substr($2,60,79) } ' >> $filelist     # Muon
    done
done




################################
# Check which files exist already
###############################
echo 'PRODUCED FILELIST LOGS'

# Copy files
echo `ls /etc/ssl/certs/CERN*pem`
cat /afs/cern.ch/user/"${USER:0:1}"/$USER/.globus/usercert.pem /etc/ssl/certs/CERN*pem >  /afs/cern.ch/user/"${USER:0:1}"/$USER/certsTemp.pem

for index in "${!PDs[@]}"
do
    filelist="filelist_$index.log"
    echo "Processing filelist: $filelist"
    cat $filelist | while read LINE
    do      
        if ! [ -e ${targetDir}${LINE:8:1000} ];
        then
            wget -r -k -4 $LINE --certificate=$certif --ca-certificate=/afs/cern.ch/user/"${USER:0:1}"/$USER/certsTemp.pem -P $targetDir 
        else
            echo "File exists already"
            echo $LINE
        fi
    done
done


#############################################
## Produce the normtag and the brilcal file
#############################################


echo 'Start with the normtag'

wget https://raw.githubusercontent.com/CMS-LUMI-POG/Normtags/master/normtag_BRIL.json
mv normtag_BRIL.json $normtagFilePath



echo 'Produce the brilcalc file'

export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH

brilcalc lumi --begin $beginRun --end $endRun -b "STABLE BEAMS" --byls -u /fb --normtag=$normtagFilePath -i /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_${beginRun}_${endRun}_Muon.json -o $brilcalcOutputFilePath

echo 'End brilcalc'


#############################################
## ZCounting step
#############################################

echo 'Start ZCounting'

mkdir -p $output_ZCounting
cd ZHarvester

./submit zmonitoring -b $beginRun -e $endRun -i /eos/cms/store/group/comm_luminosity/ZCounting/2023/DQMFiles/cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2023/ --byLsCSV $brilcalcOutputFilePath -o $output_ZCounting --etaMin $etaMin --etaCut $etaMax 


echo 'End ZCounting'
