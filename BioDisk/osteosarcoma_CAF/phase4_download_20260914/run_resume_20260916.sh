#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n run_resume_20260916.sh
# RUN_COMMAND: bash run_resume_20260916.sh
# EXPECTED_SUCCESS_MARKER: PAUSED_AFTER_DOWNLOAD_VALIDATION
set -uo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
exec 9>resume_20260916.lock
flock -n 9 || exit 73
date -Is > resume_20260916_started_at.txt
/usr/bin/time -v -o resource_usage_resume_20260916.txt /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u download_validate_resume_20260916.py --root "$PWD"
rc=$?
printf '%s\n' "$rc" > resume_20260916.exit_code
if [ "$rc" -eq 0 ] && grep -qx NON_PDX_DOWNLOAD_VALIDATION_COMPLETE COMPLETE; then
 date -Is > resume_20260916_finished_at.txt
 printf '%s\n' PAUSED_AFTER_DOWNLOAD_VALIDATION > PAUSED_AFTER_DOWNLOAD_VALIDATION
fi
exit "$rc"
