#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_preprocess_20260918
exec 8>logs/review.lock
flock -n 8 || exit 73
trap 'rc=$?; printf "%s\n" "$rc" > logs/review.exit_code' EXIT
tail --pid="$(cat logs/pipeline.pid)" -f /dev/null
[[ "$(cat logs/pipeline.exit_code)" == 0 ]]
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$PY" -m py_compile code/review_geo.py
/usr/bin/time -v -o logs/review.resources.txt "$PY" -u code/review_geo.py > logs/review.log 2>&1
