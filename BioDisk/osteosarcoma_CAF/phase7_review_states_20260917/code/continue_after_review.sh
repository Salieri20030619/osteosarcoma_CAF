#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
exec 8>logs/continuation.lock
flock -n 8 || exit 73
echo $$ > logs/continuation.pid
printf 'WAITING_FOR_REFERENCE_REVIEW %s\n' "$(date -Is)"
while kill -0 2053210 2>/dev/null; do sleep 900; done
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
# The independent identity audit is required, not merely the end of the R process.
for attempt in 1 2 3 4; do
  if grep -q REVIEW_AUDIT_READY_FOR_DAY5_6_ANALYSIS logs/report_and_audit.log 2>/dev/null && [[ -f review/final_identity_audit.json ]]; then break; fi
  sleep 30
done
if ! grep -q REVIEW_AUDIT_READY_FOR_DAY5_6_ANALYSIS logs/report_and_audit.log 2>/dev/null; then
  printf 'BLOCKED_REVIEW_AUDIT_NOT_READY %s\n' "$(date -Is)"
  exit 74
fi
if [[ -f analysis/DAY5_6_COMPUTE_COMPLETE ]]; then exit 0; fi
printf 'STARTING_DAY5_6 %s\n' "$(date -Is)"
bash code/run_day5_6.sh
printf 'COMPUTATION_COMPLETE_REQUIRES_SCIENTIFIC_REVIEW_AND_REPORT %s\n' "$(date -Is)"
