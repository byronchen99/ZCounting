#!/bin/bash

runNum=$1

workdir=$2/src/ZCounting/ZHarvester

TOP="$PWD"

echo "Prepare environment"

cd $2
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd src
eval `cmsenv`
eval `scramv1 runtime -sh`
cd $TOP

cp -r ${workdir}/Utils $TOP
cp ${workdir}/ZCounting.py $TOP

echo "run ZCounting.py"

python3 ZCounting.py -b $runNum -e $(($runNum + 1)) --input $3 --byLsCSV $4 -o $5 --mcCorrections $6 --sigModels $7 --bkgModels $8 --ptCut $9 --etaCut ${10} --mass ${11} ${12} ${13} --LumiPerMeasurement ${14} --mode ${15}
