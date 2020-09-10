#!/bin/bash

set -e

if [ "$#" -eq 1 ]; then

  if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then

    printf "\n"
    printf "%s\n" " crab3_submit_all.sh -- script to submit crab3 jobs (one per DAS dataset)"
    printf "\n"
    printf "%s\n" "  * 1 command-line argument required: path to dCache directory to host crab3 output files (crab3 parameter: \"config.Data.outLFNDirBase\")"
    printf "\n"
    exit
  fi
fi

if [ "$#" -ne 1 ]; then

  printf "\n%s\n" " >>> ERROR -- one and only one command-line argument required:"
  printf "%s\n\n" "           [1] path to dCache directory to host crab3 output files (crab3 parameter: \"config.Data.outLFNDirBase\")"
  exit
fi

TIER2_DIR="$1"

if [[ ${TIER2_DIR} != /pnfs/desy.de/cms/tier2/store/* ]]; then

  printf "\n%s\n\n" " >>> ERROR -- invalid path to dCache directory on DESY Tier-2, does not start with \"/pnfs/desy.de/cms/tier2/store/\": ${TIER2_DIR}"
  exit
fi

TIER2_DIR=${TIER2_DIR#"/pnfs/desy.de/cms/tier2"}

### --------

submit(){

 python -B \
      ${CMSSW_BASE}/src/TopAnalysis/Configuration/analysis/ttH/test/crab3/crab3_makecfg.py \
   -c ${CMSSW_BASE}/src/TopAnalysis/Configuration/analysis/common/test/ntuple_cfg.py \
   --datasets-json $(dirname $(readlink -e $0))/datasets.json \
   --tier2-dir "${TIER2_DIR}" \
   --submit \
   $@
}

### DATA ---
LUMI_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
#"${CMSSW_BASE}"/src/TopAnalysis/Configuration/python/Data/Run2017/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt

#submit  -j "${LUMI_JSON}" --DATA -s   smu_run2017B
#submit  -j "${LUMI_JSON}" --DATA -s   smu_run2017C
#submit  -j "${LUMI_JSON}" --DATA -s   smu_run2017D
#submit  -j "${LUMI_JSON}" --DATA -s   smu_run2017E
#submit  -j "${LUMI_JSON}" --DATA -s   smu_run2017F
#
#submit  -j "${LUMI_JSON}" --DATA -s    se_run2017B
#submit  -j "${LUMI_JSON}" --DATA -s    se_run2017C
#submit  -j "${LUMI_JSON}" --DATA -s    se_run2017D
#submit  -j "${LUMI_JSON}" --DATA -s    se_run2017E
#submit  -j "${LUMI_JSON}" --DATA -s    se_run2017F
#
#submit  -j "${LUMI_JSON}" --DATA -s    ee_run2017B
#submit  -j "${LUMI_JSON}" --DATA -s    ee_run2017C
#submit  -j "${LUMI_JSON}" --DATA -s    ee_run2017D
#submit  -j "${LUMI_JSON}" --DATA -s    ee_run2017E
#submit  -j "${LUMI_JSON}" --DATA -s    ee_run2017F
#
#submit  -j "${LUMI_JSON}" --DATA -s  mumu_run2017B
#submit  -j "${LUMI_JSON}" --DATA -s  mumu_run2017C
#submit  -j "${LUMI_JSON}" --DATA -s  mumu_run2017D
#submit  -j "${LUMI_JSON}" --DATA -s  mumu_run2017E
#submit  -j "${LUMI_JSON}" --DATA -s  mumu_run2017F
#
#submit  -j "${LUMI_JSON}" --DATA -s   emu_run2017B
#submit  -j "${LUMI_JSON}" --DATA -s   emu_run2017C
#submit  -j "${LUMI_JSON}" --DATA -s   emu_run2017D
#submit  -j "${LUMI_JSON}" --DATA -s   emu_run2017E
#submit  -j "${LUMI_JSON}" --DATA -s   emu_run2017F
#### --------
#
#### MC -----
#submit -s dy1050
#submit -s dy1050_ext1
#submit -s dy50_ht0070to0100
#submit -s dy50_ht0100to0200
#submit -s dy50_ht0100to0200_ext1
#submit -s dy50_ht0200to0400
#submit -s dy50_ht0200to0400_ext1
#submit -s dy50_ht0400to0600
#submit -s dy50_ht0400to0600_ext1
#submit -s dy50_ht0600to0800
#submit -s dy50_ht0800to1200
#submit -s dy50_ht1200to2500
#submit -s dy50_ht2500toINFT
#submit -s dy50inf_amcatnlofxfx
#submit -s dy50inf_amcatnlofxfx_ext1
#
#submit -s singleantitop_t
#submit -s singleantitop_tw
#submit -s singleantitop_tw_PSweights
#submit -s singletop_s
#submit -s singletop_s_PSweights
#submit -s singletop_t
#submit -s singletop_tw
#submit -s singletop_tw_PSweights
#
#submit -s tHW
#submit -s tHq
submit -s ttbarH125tobbbar
#submit -s ttbarH125tobbbar_2L
#submit -s ttbarH125tononbbbar
#
#submit -s ttbarWjetstolnu
#submit -s ttbarWjetstolnu_PSweights
#submit -s ttbarWjetstoqq
#submit -s ttbarZtollnunu
#submit -s ttbarZtoqq
#
#submit -s ttbarbg_fromDilepton
#submit -s ttbarbg_fromDilepton_PSweights
#submit -s ttbarbg_fromDilepton_matchdown
#submit -s ttbarbg_fromDilepton_matchup
#submit -s ttbarbg_fromDilepton_uetunedown
#submit -s ttbarbg_fromDilepton_uetuneup
#
#submit -s ttbarbg_fromHadronic
#submit -s ttbarbg_fromHadronic_PSweights
#submit -s ttbarbg_fromHadronic_matchdown
#submit -s ttbarbg_fromHadronic_matchup
#submit -s ttbarbg_fromHadronic_uetunedown
#submit -s ttbarbg_fromHadronic_uetuneup
#
#submit -s ttbarbg_fromLjets
#submit -s ttbarbg_fromLjets_PSweights
#submit -s ttbarbg_fromLjets_matchdown
#submit -s ttbarbg_fromLjets_matchup
#submit -s ttbarbg_fromLjets_uetunedown
#submit -s ttbarbg_fromLjets_uetuneup
#
#submit -s ttbarsignalplustau_fromDilepton
#submit -s ttbarsignalplustau_fromDilepton_PSweights
#submit -s ttbarsignalplustau_fromDilepton_matchdown
#submit -s ttbarsignalplustau_fromDilepton_matchup
#submit -s ttbarsignalplustau_fromDilepton_uetunedown
#submit -s ttbarsignalplustau_fromDilepton_uetuneup
#
#submit -s ttgjets
#submit -s ttgjets_ext1
#
#submit -s wtolnu
#submit -s wtolnu_ext1
#
#submit -s wwtoall
#submit -s wztoall
#submit -s zztoall
#### --------

unset -v TIER2_DIR
unset -v LUMI_JSON

unset -f submit
