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

cp -r ${workdir}/Utils $TOP
cp -r ${workdir}/python $TOP

echo "run ZCounting.py"
echo "./ZCounting.py -b ${runNum} -e ${runNumEnd} --input ${3} --byLsCSV ${4} -o ${5} --mcCorrections ${6} --sigModel ${7} --bkgModel ${8} --ptCut ${9} --etaCut ${10} --mass ${11} ${12} ${13} --LumiPerMeasurement ${14}"

if [[ "${15}" == "DQM" ]]; then
    cp ${workdir}/ZCounting_DQM.py $TOP    
    ./ZCounting_DQM.py -b $runNum -e ${runNumEnd} --input $3 --byLsCSV $4 -o $5 --mcCorrections $6 --sigModel $7 --bkgModel $8 --ptCut $9 --etaCut ${10} --mass ${11} ${12} ${13} --LumiPerMeasurement ${14} 
else 
    cp ${workdir}/ZCounting.py $TOP
    ./ZCounting.py -b $runNum -e ${runNumEnd} --input $3 --byLsCSV $4 -o $5 --mcCorrections $6 --sigModel $7 --bkgModel $8 --ptCut $9 --etaCut ${10} --mass ${11} ${12} ${13} --LumiPerMeasurement ${14}
fi

echo "Job done!"
