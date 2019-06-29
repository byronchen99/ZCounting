#! /bin/bash
##############
#Set Variables
#############
pass=''
cmsswRel=/afs/cern.ch/work/j/jsalfeld/Princeton/forZCounting2018/CMSSW_9_4_4/src
PD='/SingleMuon/'
COREDIR='https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2018/'
PDDIR=$COREDIR$PD
certif=/tmp/THE_CERT
targetDir=/eos/cms/store/group/comm_luminosity/ZCounting/DQMFiles2018/

##################
#Get List Of Files
#################
cd $cmsswRel
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
echo $pass | voms-proxy-init -voms cms -rfc
cd /afs/cern.ch/user/j/jsalfeld
rm -rf filelistXYZ.log
for rundir in `curl -k --cert $certif -X GET $PDDIR | awk 's=index($0,"/000") { print substr($0,s+1,10)}'`
#for rundir in `curl -k --cert $certif -X GET $PDDIR | awk ' { print substr($2,s+2,10)}'`
do
    echo ${rundir}
    curl -k --cert $certif -X GET $PDDIR$rundir | awk  -v COREDIR="$COREDIR" -F '<tr><td>' '{ print COREDIR substr($2,60,89) } ' >> filelistXYZ.log
done

################################
#Check which files exist already
###############################

echo 'PRODUCED FILELIST LOG'
##############
#Copy Files
#############
echo `ls /etc/ssl/certs/CERN*pem`
cat /afs/cern.ch/user/j/jsalfeld/.globus/usercert.pem /etc/ssl/certs/CERN*pem >  /afs/cern.ch/user/j/jsalfeld/certsTemp.pem
cat filelistXYZ.log | while read LINE
do		
	if ! [ -e ${targetDir}${LINE:8:1000} ];
	then
	wget -r -k $LINE --certificate=$certif --ca-certificate=/afs/cern.ch/user/j/jsalfeld/certsTemp.pem -P $targetDir 
	else
	echo "file exists already"
	echo $LINE
	fi
done

