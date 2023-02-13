#!/bin/bash

runNum=$1
runNumEnd=$(($1 + 1))

workdir=$2/src/ZCounting/ZHarvester

TOP="$PWD"

echo "Prepare environment"

cd $2
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd src
export HOME=pwd
export SCRAM_ARCH=el8_amd64_gcc10
eval `cmsenv`
eval `scramv1 runtime -sh`
cd $TOP

cp -r ${workdir}/Utils $TOP
cp ${workdir}/ZCounting.py $TOP

ls -lah $TOP

echo "run ZCounting.py"
echo "./ZCounting.py -b ${runNum} -e ${runNumEnd} --input ${3} --byLsCSV ${4} -o ${5} --mcCorrections ${6} --sigModel ${7} --bkgModel ${8} --ptCut ${9} --etaCut ${10} --mass ${11} ${12} ${13} --LumiPerMeasurement ${14} --mode ${15}"

./ZCounting.py -b $runNum -e ${runNumEnd} --input $3 --byLsCSV $4 -o $5 --mcCorrections $6 --sigModel $7 --bkgModel $8 --ptCut $9 --etaCut ${10} --mass ${11} ${12} ${13} --LumiPerMeasurement ${14} --mode ${15}

echo "Done running ZCounting.py"
