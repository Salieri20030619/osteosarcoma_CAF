#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n code/wait_cohort.sh
# RUN_COMMAND: bash code/wait_cohort.sh
# EXPECTED_SUCCESS_MARKER: PHASE6_COHORT_AND_INTEGRATION_FINISHED
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase6_nonpdx_cnv_20260916
pid=$(cat logs/cohort.pid)
while kill -0 "$pid" 2>/dev/null; do
 if [[ -r /proc/$pid/stat ]] && [[ $(awk '{print $3}' /proc/$pid/stat) == Z ]]; then break; fi
 sleep 900
done
cat logs/cohort.exit_code
cat STATUS.json
tail -n 12 logs/cohort.log
[[ $(cat logs/cohort.exit_code) == 0 ]]
