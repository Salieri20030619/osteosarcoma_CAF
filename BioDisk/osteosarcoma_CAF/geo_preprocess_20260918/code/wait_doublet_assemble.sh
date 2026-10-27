#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_preprocess_20260918
trap 'rc=$?; printf "%s\n" "$rc" > logs/accepted.exit_code' EXIT
tail --pid="$(cat logs/doublet.pid)" -f /dev/null
[[ "$(cat logs/doublet.exit_code)" == 0 ]]
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u code/assemble_accepted.py > logs/accepted.log 2>&1
