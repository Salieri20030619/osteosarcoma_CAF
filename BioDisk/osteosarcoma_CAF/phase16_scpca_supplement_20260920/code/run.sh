#!/usr/bin/env bash
set -uo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase16_scpca_supplement_20260920
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u code/preprocess.py > logs/preprocess.log 2>&1
rc=$?
printf '%s\n' "$rc" > logs/preprocess.exit_code
exit "$rc"
