#!/usr/bin/env bash
set -u
stage=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase11_state_reproducibility_20260920
py=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
cd "$stage" || exit 99
cp logs/analysis.log logs/analysis_initial.log
cp STATUS.json STATUS.initial.json
"$py" -m py_compile analyze.py stress_sensitivity.py || exit 98
"$py" -u stress_sensitivity.py > logs/stress_sensitivity.log 2>&1
result=$?
printf '%s\n' "$result" > stress_exit_code
if [ "$result" -ne 0 ]; then exit "$result"; fi
"$py" -u analyze.py > logs/analysis.log 2>&1
result=$?
printf '%s\n' "$result" > exit_code
exit "$result"
