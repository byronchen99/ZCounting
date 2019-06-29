#!/bin/bash
cmsswRel=/afs/cern.ch/work/j/jsalfeld/Princeton/forZCounting2018/CMSSW_9_4_4/src
cd $cmsswRel
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd ZCounting
python ZCounting.py --beginRun 315259 --parametrizeType 2
