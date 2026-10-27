#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_preprocess_20260918
exec 7>logs/doublet.lock
flock -n 7 || exit 73
trap 'rc=$?; printf "%s\n" "$rc" > logs/doublet.exit_code' EXIT
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$PY" -m py_compile code/export_doublet.py
"$PY" -u code/export_doublet.py > logs/doublet_export.log 2>&1
Rscript --vanilla code/run_scdb.R > logs/scdblfinder.log 2>&1
