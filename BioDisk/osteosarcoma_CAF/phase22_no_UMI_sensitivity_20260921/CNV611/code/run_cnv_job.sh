#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase22_no_UMI_sensitivity_20260921/CNV611
METHOD=${1:?}; SAMPLE=${2:?}
[[ "$METHOD" == infercnv || "$METHOD" == copykat ]]
[[ "$SAMPLE" =~ ^SCPCS[0-9]+$ ]]
[[ -f "$ROOT/cnv_inputs/$SAMPLE/INPUT_COMPLETE" ]]
OUT="$ROOT/$METHOD/$SAMPLE"
mkdir -p "$OUT"
exec 9>"$OUT/job.lock"
flock -n 9 || exit 73
trap 'rc=$?; printf "%s\n" "$rc" > "$OUT/job.exit_code"' EXIT
if [[ -f "$OUT/COMPLETE" ]]; then exit 0; fi
ulimit -v 67108864
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript --vanilla "$ROOT/code/run_cnv.R" "$METHOD" "$SAMPLE"
[[ -f "$OUT/COMPLETE" ]]
