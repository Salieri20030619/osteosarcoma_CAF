#!/usr/bin/env bash
set -u
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase12_study_holdout_20260920 || exit 99
mkdir -p logs
py=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$py" -m py_compile run_validation.py || exit 98
"$py" -u run_validation.py > logs/validation.log 2>&1
result=$?
printf '%s\n' "$result" > exit_code
exit "$result"
