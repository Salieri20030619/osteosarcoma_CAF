#!/usr/bin/env bash
set -u
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase13_continuous_programs_20260920 || exit 99
mkdir -p logs
py=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$py" -m py_compile continuous.py || exit 98
"$py" -u continuous.py > logs/continuous.log 2>&1
result=$?
printf '%s\n' "$result" > exit_code
exit "$result"
