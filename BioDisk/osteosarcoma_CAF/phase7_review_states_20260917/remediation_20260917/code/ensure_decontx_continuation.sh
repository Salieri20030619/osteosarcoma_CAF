#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
exec 7>remediation_20260917/logs/export_continuation.lock
flock -n 7 || exit 73
while kill -0 2418610 2>/dev/null; do sleep 900; done
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u remediation_20260917/code/export_decontx.py > remediation_20260917/logs/export_resume.log 2>&1
if [[ ! -f remediation_20260917/logs/DECONTX_RUNS_EXITED ]]; then
 set +e
 bash remediation_20260917/code/run_decontx.sh
 rc=$?
 if [[ $rc -ne 0 && $rc -ne 73 ]]; then exit "$rc"; fi
fi
