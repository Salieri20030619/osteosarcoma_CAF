#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n code/observe_completion.sh
# RUN_COMMAND: bash code/observe_completion.sh
# EXPECTED_SUCCESS_MARKER: PHASE6_COMPLETION_STATE_RECORDED
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase6_nonpdx_cnv_20260916
pid=$(cat logs/cohort.pid)
while kill -0 "$pid" 2>/dev/null; do
 if [[ -r /proc/$pid/stat ]] && [[ $(awk '{print $3}' /proc/$pid/stat) == Z ]]; then break; fi
 sleep 900
done
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python code/record_completion.py
