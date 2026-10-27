#!/usr/bin/env bash
set -uo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_identity_20260919
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
exec 8>logs/pipeline.lock
flock -n 8 || exit 73
echo $$ > logs/pipeline.pid
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u code/prepare_and_run.py > logs/pipeline.log 2>&1
rc=$?
printf '%s\n' "$rc" > logs/pipeline.exit_code
exit "$rc"
