#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n code/finish_qc_review.sh
# RUN_COMMAND: bash code/finish_qc_review.sh
# EXPECTED_SUCCESS_MARKER: PHASE5_REVIEW_AND_AUDIT_COMPLETE
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase5_nonpdx_qc_20260916
pid=$(cat logs/all.pid)
while kill -0 "$pid" 2>/dev/null; do
 if [[ -r /proc/$pid/stat ]] && [[ $(awk '{print $3}' /proc/$pid/stat) == Z ]]; then break; fi
 sleep 900
done
[[ $(cat logs/all.exit_code) == 0 ]]
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
/usr/bin/time -v -o logs/review.resources.txt "$PY" -u code/review_lineages.py > logs/review.log 2>&1
/usr/bin/time -v -o logs/audit.resources.txt "$PY" -u code/validate_stage.py > logs/audit.log 2>&1
echo PHASE5_REVIEW_AND_AUDIT_COMPLETE
