#!/usr/bin/env bash
set -euo pipefail
OUT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase90_patient_validation_20261015/cnv/infercnv/"$1"
if [[ -f /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase90_patient_validation_20261015/DEFER_NEW_HMM.json && ! -f "$OUT/COMPLETE" ]]; then
  echo HMM_DEFERRED_AFTER_REFERENCE_VALIDATION_FAILURE >&2
  exit 75
fi
mkdir -p "$OUT"
exec 9>"$OUT/phase90_execution.lock"
flock -x 9
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
ulimit -v 100663296
exec Rscript --vanilla /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase90_patient_validation_20261015/code/native_cnv_v2_cleanrefs.R infercnv "$1"
