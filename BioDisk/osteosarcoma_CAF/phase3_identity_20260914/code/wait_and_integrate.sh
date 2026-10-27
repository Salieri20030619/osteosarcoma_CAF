#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n wait_and_integrate.sh
# RUN_COMMAND: bash wait_and_integrate.sh
# EXPECTED_SUCCESS_MARKER: IDENTITY_PIPELINE_COMPLETE
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
echo $$ > "$ROOT/logs/wait_and_integrate.pid"
files=("$ROOT/copykat/S0058/job.pid" "$ROOT/logs/cohort_infercnv_0.pid" "$ROOT/logs/cohort_copykat_1.pid" "$ROOT/logs/cohort_copykat_2.pid")
for file in "${files[@]}"; do test -f "$file"; done
while :; do
  running=0
  for file in "${files[@]}"; do
    pid=$(cat "$file");[[ "$pid" =~ ^[0-9]+$ ]]
    if kill -0 "$pid" 2>/dev/null; then
      state=$(ps -p "$pid" -o stat= || true)
      if [[ -n "$state" && "$state" != Z* ]]; then running=1; fi
    fi
  done
  if (( running == 0 )); then break; fi
  sleep 900
done
failed=0
for method in infercnv copykat; do
  for sample in S0058 S0059 S0217 S0218 SC069 SC072 SC073; do
    if [[ ! -f "$ROOT/$method/$sample/COMPLETE" ]]; then echo "CNV_NOT_COMPLETE $method $sample";failed=1;fi
  done
done
if (( failed != 0 )); then echo CNV_COMPLETION_CHECK_FAILED;exit 1;fi
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$PY" -m py_compile "$ROOT/code/integrate_identity.py"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MPLBACKEND=Agg
/usr/bin/time -v -o "$ROOT/logs/identity_resource_usage.txt" "$PY" "$ROOT/code/integrate_identity.py" --root /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/OS_CAF_phase2 --phase3 "$ROOT"
echo IDENTITY_PIPELINE_COMPLETE
