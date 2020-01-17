#!/bin/bash

runNum=$1

workdir=$2/src/ZCounting/ZHarvester

TOP="$PWD"

cd $2
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd src
eval `cmsenv`
eval `scramv1 runtime -sh`
cd $TOP

cp -r ${workdir}/Utils $TOP
cp -r ${workdir}/.rootlogon.C $TOP
cp ${workdir}/calculateDataEfficiency.C $TOP
cp ${workdir}/calculateZEfficiency.C $TOP
cp ${workdir}/ZCounting.py $TOP

python ZCounting.py -b $runNum -e $(($runNum + 1)) --dirDQM $3 --byLsCSV $4 -o $5 --mcCorrections $6 --sigTemplates $7 --ptCut $8
