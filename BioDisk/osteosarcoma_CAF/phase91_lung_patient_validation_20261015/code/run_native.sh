#!/usr/bin/env bash
set -euo pipefail
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
ulimit -v 100663296
exec Rscript --vanilla /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase91_lung_patient_validation_20261015/code/native_cnv.R infercnv "$1"
