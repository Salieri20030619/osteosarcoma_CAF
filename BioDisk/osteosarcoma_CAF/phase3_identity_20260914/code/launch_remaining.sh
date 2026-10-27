#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n launch_remaining.sh
# RUN_COMMAND: bash launch_remaining.sh
# EXPECTED_SUCCESS_MARKER: CNV_REMAINING_LAUNCHED
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
bash -n "$ROOT/code/run_cnv_cohort.sh"
for pair in infercnv:0 copykat:1 copykat:2; do
  method=${pair%:*}; shard=${pair#*:}
  pidfile="$ROOT/logs/cohort_${method}_${shard}.pid"
  if [[ -f "$pidfile" ]] && kill -0 "$(cat "$pidfile")" 2>/dev/null; then echo "ALREADY_RUNNING $pair";continue;fi
  nohup bash "$ROOT/code/run_cnv_cohort.sh" "$method" "$shard" > "$ROOT/logs/cohort_${method}_${shard}.log" 2>&1 < /dev/null &
  printf '%s\t%s\n' "$pair" "$!" >> "$ROOT/logs/cohort_launch.tsv"
done
echo CNV_REMAINING_LAUNCHED
