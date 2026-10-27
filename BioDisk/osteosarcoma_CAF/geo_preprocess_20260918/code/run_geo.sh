#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_preprocess_20260918
exec 9>logs/pipeline.lock
flock -n 9 || exit 73
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4 MPLBACKEND=Agg
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
trap 'rc=$?; printf "%s\n" "$rc" > logs/pipeline.exit_code; date -u +%FT%TZ > logs/finished_utc.txt' EXIT
date -u +%FT%TZ > logs/started_utc.txt
"$PY" -m py_compile code/prepare_geo.py
/usr/bin/time -v -o logs/prepare.resources.txt "$PY" -u code/prepare_geo.py > logs/prepare.log 2>&1
"$PY" -m py_compile code/geo_qc.py
/usr/bin/time -v -o logs/qc.resources.txt "$PY" -u code/geo_qc.py --manifest manifest.tsv --raw inputs --out results > logs/qc.log 2>&1
grep -qx GEO_QC_CLUSTERING_COMPLETE results/SUCCESS.txt
printf '%s\n' GEO_PREPROCESS_COMPUTATION_COMPLETE > COMPUTATION_COMPLETE
