#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase88_expansion_identity_20261014
METHOD=${1:?}; SAMPLE=${2:?}
[[ "$METHOD" == infercnv || "$METHOD" == copykat || "$METHOD" == decontx ]]
[[ "$SAMPLE" =~ ^[A-Za-z0-9_]+$ ]]
OUT="$ROOT/$METHOD/$SAMPLE"
mkdir -p "$OUT"
exec 9>"$OUT/job.lock"
flock -n 9 || exit 73
echo $$ > "$OUT/job.pid"
trap 'rc=$?; printf "%s\n" "$rc" > "$OUT/job.exit_code"' EXIT
ulimit -v 67108864
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
if [[ "$METHOD" == decontx ]]; then
 /usr/bin/time -v -o "$OUT/resource_usage.txt" /mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla "$ROOT/code/run_decontx.R" "$SAMPLE"
else
 source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
 conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
 /usr/bin/time -v -o "$OUT/resource_usage.txt" Rscript --vanilla "$ROOT/code/run_cnv.R" "$METHOD" "$SAMPLE"
fi
[[ -f "$OUT/COMPLETE" ]]
