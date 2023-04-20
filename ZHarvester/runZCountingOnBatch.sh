#!/bin/bash

runNum=$1
runNumEnd=$(($1 + 1))

workdir=$2/src/ZCounting/ZHarvester

TOP="$PWD"

cd $2
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd src
eval `cmsenv`
eval `scramv1 runtime -sh`
cd $TOP

cp -r ${workdir}/roofit $TOP
cp -r ${workdir}/common $TOP


if [[ "${16}" == "DQM" ]]; then
    echo "./ZCounting_DQM.py -b ${runNum} -e ${runNumEnd} --input ${3} --byLsCSV ${4} -o ${5} --mcCorrections ${6} --sigModel ${7} --bkgModel ${8} --ptCut ${9} --etaMin ${10} --etaCut ${11} --mass ${12} ${13} ${14} --LumiPerMeasurement ${15}"

    cp ${workdir}/ZCounting_DQM.py $TOP    
    ./ZCounting_DQM.py -b $runNum -e ${runNumEnd} --input $3 --byLsCSV $4 -o $5 --mcCorrections $6 --sigModel $7 --bkgModel $8 --ptCut $9 --etaMin ${10} --etaCut ${11} --mass ${12} ${13} ${14} --LumiPerMeasurement ${15} 
else 
    echo "./ZCounting.py -b ${runNum} -e ${runNumEnd} --input ${3} --byLsCSV ${4} -o ${5} --mcCorrections ${6} --sigModel ${7} --bkgModel ${8} --ptCut ${9} --etaMin ${10} --etaCut ${11} --mass ${12} ${13} ${14} --LumiPerMeasurement ${15}"
    cp ${workdir}/ZCounting.py $TOP
    ./ZCounting.py -b $runNum -e ${runNumEnd} --input $3 --byLsCSV $4 -o $5 --mcCorrections $6 --sigModel $7 --bkgModel $8 --ptCut $9 --etaMin ${10} --etaCut ${11} --mass ${12} ${13} ${14} --LumiPerMeasurement ${15}
fi
