#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n wait_resume_20260916.sh
# RUN_COMMAND: bash wait_resume_20260916.sh
# EXPECTED_SUCCESS_MARKER: PAUSED_AFTER_DOWNLOAD_VALIDATION
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
job_pid=$(cat resume_20260916.pid)
while kill -0 "$job_pid" 2>/dev/null; do
 if [[ -r /proc/$job_pid/stat ]] && [[ $(awk '{print $3}' /proc/$job_pid/stat) == Z ]]; then break; fi
 sleep 900
done
cat resume_20260916.exit_code
tail -n 30 download_resume_20260916.log
test "$(cat resume_20260916.exit_code)" = 0
grep -qx NON_PDX_DOWNLOAD_VALIDATION_COMPLETE COMPLETE
cat PAUSED_AFTER_DOWNLOAD_VALIDATION
