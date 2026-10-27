#!/usr/bin/env bash
set -euo pipefail
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript --vanilla /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918/code/preflight.R
