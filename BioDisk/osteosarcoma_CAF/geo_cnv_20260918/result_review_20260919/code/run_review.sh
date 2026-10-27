#!/usr/bin/env bash
set -euo pipefail
P=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918
O="$P/result_review_20260919"
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
"$PY" "$O/code/review.py" > "$O/logs/review.log" 2>&1
"$PY" "$O/code/reference_review.py" > "$O/logs/reference_review.log" 2>&1
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
Rscript --vanilla "$O/code/export_profiles.R" "$P" > "$O/logs/profile_export.log" 2>&1
printf '0\n' > "$O/logs/profile_export.exit_code"
"$PY" "$O/code/profile_review.py" > "$O/logs/profile_review.log" 2>&1
printf '0\n' > "$O/logs/profile_review.exit_code"
# Reproduction stops before finalization; inspect new figures, then run finalize_review.py.
printf 'REVIEW_TABLES_AND_FIGURES_REGENERATED_VISUAL_QA_REQUIRED\n'
