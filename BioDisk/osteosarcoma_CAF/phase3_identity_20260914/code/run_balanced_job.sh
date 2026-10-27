#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n run_balanced_job.sh
# RUN_COMMAND: nohup bash run_balanced_job.sh <S0058|SC072> > <log> 2>&1 < /dev/null &
# EXPECTED_SUCCESS_MARKER: BALANCED_JOB_COMPLETE
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
SAMPLE=${1:?sample}
[[ "$SAMPLE" =~ ^(S0058|SC072)$ ]]
OUT="$ROOT/copykat_balanced/$SAMPLE"
mkdir -p "$OUT"
test ! -f "$OUT/COMPLETE" || { echo BALANCED_JOB_CACHED;exit 0; }
mkdir "$OUT/job.lock"
echo $$ > "$OUT/job.pid"
date --iso-8601=seconds > "$OUT/job.started"
trap 'rc=$?; printf "%s\n" "$rc" > "$OUT/job.exit_code"; rmdir "$OUT/job.lock"' EXIT
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' "$ROOT/code/run_copykat_balanced.R"
/usr/bin/time -v -o "$OUT/resource_usage.txt" Rscript --vanilla "$ROOT/code/run_copykat_balanced.R" copykat "$SAMPLE"
test -f "$OUT/COMPLETE"
echo BALANCED_JOB_COMPLETE
