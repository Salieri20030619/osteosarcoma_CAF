#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n code/run_cnv_job.sh
# RUN_COMMAND: bash code/run_cnv_job.sh <infercnv|copykat> <sample>
# EXPECTED_SUCCESS_MARKER: CNV_JOB_COMPLETE
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
METHOD=scevan; SAMPLE=${1:?}

[[ "$SAMPLE" =~ ^[A-Za-z0-9_.]+$ ]]
[[ -f "$ROOT/inputs/$SAMPLE/INPUT_COMPLETE" ]]
OUT="$ROOT/$METHOD/$SAMPLE"
mkdir -p "$OUT"
exec 9>"$OUT/job.lock"
flock -n 9 || exit 73
[[ ! -f "$OUT/COMPLETE" ]] || { echo CNV_JOB_CACHED; exit 0; }
echo $$ > "$OUT/job.pid"
date -Is > "$OUT/job.started"
trap 'rc=$?; printf "%s\n" "$rc" > "$OUT/job.exit_code"' EXIT
# Resource gate is checked in the shell, without repeated model polling.
while [[ $(awk '/MemAvailable/{print $2}' /proc/meminfo) -lt 67108864 ]]; do sleep 900; done
ulimit -v 67108864
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
Rscript --vanilla "$ROOT/code/check_R_syntax.R" "$ROOT/code/run_scevan.R"
/usr/bin/time -v -o "$OUT/resource_usage.txt" Rscript --vanilla "$ROOT/code/run_scevan.R" "$SAMPLE"
[[ -f "$OUT/COMPLETE" ]]
echo CNV_JOB_COMPLETE
