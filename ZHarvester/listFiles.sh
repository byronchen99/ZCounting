#!/bin/bash
PPHRASE="YOUR_PPHRASE"

PD='SingleMuon/'
COREDIR='https://cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2017/'
PDDIR=$COREDIR$PD

for rundir in `curl -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem  --pass "YOUR_PPHRASE" -X GET $PDDIR | awk 's=index($0,"/000") { print substr($0,s+1,10)}'` 
do
    curl -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem --pass $PPHASE -X GET $PDDIR$rundir | awk  -v COREDIR="$COREDIR" -F '<tr><td>' '{ print COREDIR substr($2,60,89) } ' >> filelist.log
done

