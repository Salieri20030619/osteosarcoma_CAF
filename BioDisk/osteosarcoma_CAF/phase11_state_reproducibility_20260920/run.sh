#!/usr/bin/env bash
set -u
stage=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase11_state_reproducibility_20260920
cd "$stage" || exit 99
mkdir -p logs
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -m py_compile analyze.py || exit 98
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u analyze.py > logs/analysis.log 2>&1
result=$?
printf '%s\n' "$result" > exit_code
exit "$result"
