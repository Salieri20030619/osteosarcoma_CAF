#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase92_copykat_validation_20261015
mkdir -p "$ROOT/copykat/$1"
exec 9>"$ROOT/copykat/$1/execution.lock"
flock -n 9
ulimit -v 100663296
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
exec Rscript --vanilla "$ROOT/code/copykat_full.R" "$ROOT" "$1"
