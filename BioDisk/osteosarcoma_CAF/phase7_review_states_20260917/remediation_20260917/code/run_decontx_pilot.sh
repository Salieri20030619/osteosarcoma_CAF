#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate R_4.5.2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
flock remediation_20260917/logs/decontx_GSE152048__BC10.lock Rscript --vanilla remediation_20260917/code/run_decontx.R GSE152048__BC10
