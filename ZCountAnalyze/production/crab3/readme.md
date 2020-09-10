
crab3 interface
===============

* **Step 0**: `source crab3_env.sh`

* **Step 0.1**: if starting a new production, `crab3_localarea.py --mask-data` (see script description below)

* **Step 0.2**: if starting a new production, `crab3_clean_cache.py` (see script description below)

* **Step 1**: prepare the **datasets JSON**, i.e. a JSON file containing the configuration of the TopAnalysis NTuples to be produced via crab3.
              This file needs to be prepared by hand.
              Examples can be found in the `datasets/` directory.

* **Step 2**: run `crab3_production_json.py` to create the **production JSON**.
              The latter contains all the information needed for submitting a set of crab3 tasks.
              Necessary inputs include the datasets JSON(s), the target Tier-2 storage area and the target output directory for the final TopAnalysis NTuples.
              This JSON file should be kept at hand as long as the corresponding tasks are running (see `crab3_monitor.py`).
              For more info, do `crab3_production_json.py -h`.

* **Step 3**: run `crab3_submit.py` to submit crab3 tasks.
              Each task is defined as an entry in the **production JSON** used as input to `crab3_submit.py`.
              For more info, do `crab3_submit.py -h`.

* **Step 4**: run `crab3_monitor.py` periodically to check the status of the specified crab3 tasks.
              Option `--resubmit`: failed jobs will be resubmitted.
              Option `--hadd`: the outputs of a fully completed task are merged into the final TopAnalysis NTuple.
              Some options, e.g. `--hadd`, require the production JSON as input.
              For more info, do `crab3_monitor.py -h`.

* **Step 5**: when re-logging in to check on the tasks, do `source crab3_env.sh` (**Step 0**), then move directly to **Step 4**.

Utilities
===============================

* **utilities/** directory: contains scripts to facilitate simple routine tasks

  - `crab3_env.sh`:
     * source this file to set up CMSSW and crab3 environments

  - `crab3_localarea.py --mask-data`:
     * WARNING: executing this script will rename some directories, double-check before executing
     * renames `data/` sub-directories in `$CMSSW_BASE` that are not required for NTuple production
     * necessary to bring size of crab3 tarball below max-size limit
     * this script is supposed to be executed once, before submitting jobs
     * if not executed, grid jobs will most likely fail due to crab3 tarball exceeding max-size limit

  - `crab3_clean_cache.py`:
     * cleans up the CRAB User File Cache
     * script taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ#User_quota_in_the_CRAB_cache

  - `crab3_all.sh`:
     * simple wrapper to apply the same crab3 command to all the crab3 folders in the current directory
     * example: `./crab3_all.sh resubmit`
