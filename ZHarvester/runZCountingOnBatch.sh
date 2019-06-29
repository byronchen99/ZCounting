#!/bin/bash

runNum=$1

workdir=$2/src/ZCounting/

TOP="$PWD"

cd $2
eval `scramv1 runtime -sh`
cd $TOP

cp -r ${workdir}/Utils $TOP
cp -r ${workdir}/.rootlogon.C $TOP
cp ${workdir}/calculateDataEfficiency.C $TOP
cp ${workdir}/calculateZEfficiency.C $TOP
cp ${workdir}/ZCounting.py $TOP

python ZCounting.py -b $runNum -e $(($runNum + 1)) -d $3 -f $4 -a $5 -x $6

