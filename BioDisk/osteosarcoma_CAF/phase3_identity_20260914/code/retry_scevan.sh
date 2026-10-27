#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n retry_scevan.sh
# RUN_COMMAND: bash retry_scevan.sh
# EXPECTED_SUCCESS_MARKER: SCEVAN_INTERFACE_RETRY_LAUNCHED
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
for s in S0058 SC072; do
  OUT="$ROOT/scevan/$s"
  test ! -d "$OUT/job.lock"
  if [[ -f "$OUT/FAILED" ]]; then mv "$OUT/FAILED" "$OUT/FAILED_interface_attempt.txt";fi
  if [[ -f "$ROOT/logs/scevan_$s.log" ]]; then mv "$ROOT/logs/scevan_$s.log" "$ROOT/logs/scevan_${s}_interface_attempt.log";fi
  nohup bash "$ROOT/code/run_scevan_job.sh" "$s" > "$ROOT/logs/scevan_$s.log" 2>&1 < /dev/null &
done
echo SCEVAN_INTERFACE_RETRY_LAUNCHED
