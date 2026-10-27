#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n code/run_day5_6.sh
# RUN_COMMAND: bash code/run_day5_6.sh
# EXPECTED_SUCCESS_MARKER: ORIGINAL_DAY5_6_COMPUTE_COMPLETE
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
exec 9>logs/day5_6.lock
flock -n 9 || exit 73
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
[[ -f review/final_identity_audit.json ]]
grep -q REVIEW_AUDIT_READY_FOR_DAY5_6_ANALYSIS logs/report_and_audit.log
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4 MPLBACKEND=Agg
echo $$ > logs/day5_6.pid
date -Is > logs/day5_6.started
trap 'rc=$?; printf "%s\n" "$rc" > logs/day5_6.exit_code' EXIT
/usr/bin/time -v -o logs/states.resources.txt "$PY" -u code/analyze_states.py > logs/states.log 2>&1
/usr/bin/time -v -o logs/mycaf.resources.txt "$PY" -u code/mycaf_effects.py > logs/mycaf.log 2>&1
/usr/bin/time -v -o logs/lineage_programs.resources.txt "$PY" -u code/compare_lineage_programs.py > logs/lineage_programs.log 2>&1
"$PY" -u code/audit_day5_6_design.py > logs/design_audit.log 2>&1
"$PY" -u code/audit_original_day1_6.py > logs/plan_audit.log 2>&1
printf '%s\n' ORIGINAL_DAY5_6_COMPUTE_COMPLETE > analysis/DAY5_6_COMPUTE_COMPLETE
date -Is > logs/day5_6.finished
