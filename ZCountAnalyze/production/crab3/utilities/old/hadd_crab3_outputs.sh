#!/bin/bash

if [ "$#" -eq 1 ]; then

  if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then

    printf "\n"
    printf "%s\n" " hadd_crab3_outputs.sh -- script to merge crab3 outputs"
    printf "\n"
    printf "%s\n" "  * 2 command-line arguments required:"
    printf "%s\n" "      [1] path to dCache directory specified in crab3 configuration file(s)"
    printf "%s\n" "      [2] path to directory for final output NTuples"
    printf "\n"
    exit
  fi
fi

if [ "$#" -ne 2 ]; then

  printf "\n%s\n" " >>> ERROR -- invalid command-line argument(s):"
  printf "%s\n"   "           [1] path to dCache directory specified in crab3 configuration file(s)"
  printf "%s\n\n" "           [2] path to directory for final output NTuples"
  exit
fi

IDIR="$1"
ODIR="$2"

### --------

if [ ! -d "${IDIR}" ]; then

  printf "\n%s\n\n" " >>> ERROR -- path to input dCache directory not found: ${IDIR}"
  exit
fi

hadd_outputs(){

  if [ "$#" -ne 4 ]; then return 1; fi;

  local EXE="${CMSSW_BASE}/src/TopAnalysis/Configuration/analysis/common/test/scripts/hadd_wrapper.py"

  if [ ! -f "${EXE}" ]; then return 1; fi;

  $EXE "$4".root "$1"/"$4"/*/*/*/*/"$4"_*.root -s "$2" -N "$3"
}

mkdir -p "${ODIR}"
cd       "${ODIR}"

### --------

hadd_outputs  "${IDIR}"   1   158  ee_run2017B
hadd_outputs  "${IDIR}"   1   257  ee_run2017C
hadd_outputs  "${IDIR}"   1    33  ee_run2017D
hadd_outputs  "${IDIR}"   1   129  ee_run2017E
hadd_outputs  "${IDIR}"   1   172  ee_run2017F

hadd_outputs  "${IDIR}"   1    35  emu_run2017B
hadd_outputs  "${IDIR}"   1    51  emu_run2017C
hadd_outputs  "${IDIR}"   1    22  emu_run2017D
hadd_outputs  "${IDIR}"   1   133  emu_run2017E
hadd_outputs  "${IDIR}"   1    68  emu_run2017F

hadd_outputs  "${IDIR}"   1    21  mumu_run2017B
hadd_outputs  "${IDIR}"   1   154  mumu_run2017C
hadd_outputs  "${IDIR}"   1    33  mumu_run2017D
hadd_outputs  "${IDIR}"   1   147  mumu_run2017E
hadd_outputs  "${IDIR}"   1   232  mumu_run2017F

hadd_outputs  "${IDIR}"   1   191  se_run2017B
hadd_outputs  "${IDIR}"   2  1182  se_run2017C
hadd_outputs  "${IDIR}"   1   337  se_run2017D
hadd_outputs  "${IDIR}"   1   572  se_run2017E
hadd_outputs  "${IDIR}"   1   574  se_run2017F

hadd_outputs  "${IDIR}"   1    46  smu_run2017B
hadd_outputs  "${IDIR}"   1   162  smu_run2017C
hadd_outputs  "${IDIR}"   1    84  smu_run2017D
hadd_outputs  "${IDIR}"   1   155  smu_run2017E
hadd_outputs  "${IDIR}"   2   798  smu_run2017F

hadd_outputs  "${IDIR}"   1    66  ttbarsignalplustau_fromDilepton
hadd_outputs  "${IDIR}"  10   891  ttbarsignalplustau_fromDilepton_PSweights

hadd_outputs  "${IDIR}"   1    27  ttbarbg_fromDilepton
hadd_outputs  "${IDIR}"   2   214  ttbarbg_fromDilepton_PSweights
hadd_outputs  "${IDIR}"   1    70  ttbarbg_fromLjets
hadd_outputs  "${IDIR}"   5   410  ttbarbg_fromLjets_PSweights

hadd_outputs  "${IDIR}"   1     4  singleantitop_t
hadd_outputs  "${IDIR}"   1    23  singleantitop_tw
hadd_outputs  "${IDIR}"   1    14  singleantitop_tw_PSweights
hadd_outputs  "${IDIR}"   1    39  singletop_s
hadd_outputs  "${IDIR}"   1    19  singletop_s_PSweights
hadd_outputs  "${IDIR}"   1    27  singletop_t
hadd_outputs  "${IDIR}"   1    21  singletop_tw
hadd_outputs  "${IDIR}"   1    53  singletop_tw_PSweights

hadd_outputs  "${IDIR}"   1    12  dy1050
hadd_outputs  "${IDIR}"   1   153  dy50inf_amcatnlofxfx
hadd_outputs  "${IDIR}"  10  1083  dy50inf_amcatnlofxfx_ext1

hadd_outputs  "${IDIR}"   1    25  wwtoall
hadd_outputs  "${IDIR}"   1     4  wztoall
hadd_outputs  "${IDIR}"   1     2  zztoall

hadd_outputs  "${IDIR}"   1    27  ttbarWjetstolnu
hadd_outputs  "${IDIR}"   1    42  ttbarWjetstolnu_PSweights
hadd_outputs  "${IDIR}"   1     5  ttbarWjetstoqq
hadd_outputs  "${IDIR}"   1   188  ttbarZtollnunu
hadd_outputs  "${IDIR}"   1    23  ttbarZtoqq
hadd_outputs  "${IDIR}"   1    13  ttgjets

hadd_outputs  "${IDIR}"   1   342  ttbarH125tobbbar
hadd_outputs  "${IDIR}"   5   101  ttbarH125tobbbar_2L
hadd_outputs  "${IDIR}"   1   264  ttbarH125tononbbbar

### --------

cd "${OLDPWD}"

unset -f hadd_outputs
unset -v IDIR ODIR
