#! /bin/bash

STORAGE=$1
PD='SingleMuon/'
COREDIR='https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2017/'
PDDIR=$COREDIR$PD

echo ">>> enter passphrase ..."
read -sp 'passphrase: ' PPHRASE 

echo ">>> scan files on cmsweb"
rm filelist.log
for rundir in `curl -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem  --pass $PPHRASE -X GET $PDDIR | awk 's=index($0,"/000") { print substr($0,s+1,10)}'` 
do
    ISTRING=$(curl -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem --pass $PPHRASE -X GET $PDDIR$rundir | awk  -v COREDIR="$COREDIR" -F '<tr><td>' '{ print COREDIR substr($2,60,index($2,".root")-60+5) } ')

    for substring in ${ISTRING}
    do
        if [[ $substring != *09Aug2019_UL2017-v1__DQMIO.root ]] ;
        then
            continue
        fi
        echo $substring >> filelist.log
    done
done

echo ">>> load files"
cat ~/.globus/usercert.pem /etc/ssl/certs/CERN*pem >  /tmp/$USER/certs.pem
certU=/tmp/x509up_u30079
cat filelist.log | while read LINE
do
wget -r -k $LINE --certificate=$certU --ca-certificate=/tmp/$USER/certs.pem -P $STORAGE 
echo $LINE
done
