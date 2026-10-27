#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n code/run_qc.sh
# RUN_COMMAND: bash code/run_qc.sh pilot OR bash code/run_qc.sh all
# EXPECTED_SUCCESS_MARKER: PHASE5_PILOTS_COMPLETE OR PHASE5_QC_CLUSTERING_COMPLETE
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase5_nonpdx_qc_20260916
mode=${1:?pilot or all}
[[ "$mode" == pilot || "$mode" == all ]]
exec 9>logs/qc.lock
flock -n 9 || exit 73
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4 MPLBACKEND=Agg
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$PY" -m py_compile code/phase5.py
trap 'rc=$?; printf "%s\n" "$rc" > logs/'"$mode"'.exit_code' EXIT
if [[ "$mode" == pilot ]]; then
 for sample in GSE299023__OS.C.21.BX.OCT GSE152048__BC17; do
  /usr/bin/time -v -o "logs/${sample}.resources.txt" "$PY" -u code/phase5.py --manifest manifest.tsv --raw ../phase4_download_20260914/matrices --out results --sample "$sample" > "logs/${sample}.log" 2>&1
 done
 printf '%s\n' PHASE5_PILOTS_COMPLETE > PILOTS_COMPLETE
else
 /usr/bin/time -v -o logs/all.resources.txt "$PY" -u code/phase5.py --manifest manifest.tsv --raw ../phase4_download_20260914/matrices --out results > logs/all.log 2>&1
 grep -qx PHASE5_QC_CLUSTERING_COMPLETE results/SUCCESS.txt
fi
