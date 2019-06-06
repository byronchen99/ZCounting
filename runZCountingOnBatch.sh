#!/bin/bash

runNum=$1

workdir=${CMSSW_BASE}/src/ZCounting/

TOP="$PWD"

cd $CMSSW_BASE
eval `scramv1 runtime -sh`
cd $TOP

cp -r ${workdir}/Utils $TOP
cp -r ${workdir}/.rootlogon.C $TOP
cp ${workdir}/calculateDataEfficiency.C $TOP
cp ${workdir}/calculateZEfficiency.C $TOP
cp ${workdir}/ZCounting.py $TOP

python ZCounting.py -b $runNum -e $(($runNum + 1)) -d $2 -f $3 -g $4 -t $5 -a $6 -x $7  

