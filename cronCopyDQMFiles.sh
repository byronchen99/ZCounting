#! /bin/bash

##############
#Set Variables
#############

pass='rainbowsix'            # grid passwd
cmsswRel='/afs/cern.ch/user/t/tmenezes/work/private/CMSSW_12_4_2/src/'   # path to ZCounting dir.
#PD='/SingleMuon/'
# After the RunID = 356426, SingleMuon -> Muon 
PD='/Muon/'
COREDIR='https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2022/'
PDDIR=$COREDIR$PD
certif=/tmp/x509up_u81013     # you get one after the voms (grid proxy)
targetDir=/eos/cms/store/group/comm_luminosity/ZCounting/2022/DQMFiles/

##################
#Get List Of Files
#################

cd $cmsswRel
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
echo $pass | voms-proxy-init -voms cms -rfc    # get the grid proxy 
cd /afs/cern.ch/user/t/tmenezes/work/private/CMSSW_12_4_2/src/ZCounting  # enter the ZCounting dir.
rm -rf filelist.log


# Create the log files with the files from the PDDIR (uses grid certs. and grid passwd)
for rundir in `curl -k --cert $certif -X GET $PDDIR | awk 's=index($0,"/000") { print substr($0,s+1,10)}'` 
do
    #curl -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem --pass $pass -X GET $PDDIR$rundir | awk  -v COREDIR="$COREDIR" -F '<tr><td>' '{ print COREDIR substr($2,60,89) } ' >> filelist.log    # SingleMuon
    curl -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem --pass $pass -X GET $PDDIR$rundir | awk  -v COREDIR="$COREDIR" -F '<tr><td>' '{ print COREDIR substr($2,60,77) } ' >> filelist.log     # Muon
done


################################
#Check which files exist already
###############################
echo 'PRODUCED FILELIST LOG'

# Copy files
echo `ls /etc/ssl/certs/CERN*pem`
cat /afs/cern.ch/user/t/tmenezes/.globus/usercert.pem /etc/ssl/certs/CERN*pem >  /afs/cern.ch/user/t/tmenezes/certsTemp.pem 
cat filelist.log | while read LINE
do		
	if ! [ -e ${targetDir}${LINE:8:1000} ];
    then
    wget -r -k -4 $LINE --certificate=$certif --ca-certificate=/afs/cern.ch/user/t/tmenezes/certsTemp.pem -P $targetDir  
    else
	echo "file exists already"
	echo $LINE
	fi
done
