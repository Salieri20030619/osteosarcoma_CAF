#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n wait_cnv_jobs.sh
# RUN_COMMAND: bash wait_cnv_jobs.sh <job.pid> [job.pid ...]
# EXPECTED_SUCCESS_MARKER: COMPLETION_WAIT_FINISHED
# A completion-only wait: process existence checked every 15 minutes, no repeated model polling.
set -euo pipefail
for file in "$@"; do
  test -f "$file"
  pid=$(cat "$file")
  [[ "$pid" =~ ^[0-9]+$ ]]
  while kill -0 "$pid" 2>/dev/null; do
    state=$(ps -p "$pid" -o stat= || true)
    [[ "$state" != Z* ]] || break
    sleep 900
  done
  dir=$(dirname "$file")
  if [[ -f "$dir/COMPLETE" ]]; then printf 'COMPLETE %s\n' "$dir"; else printf 'NOT_COMPLETE %s\n' "$dir"; fi
done
echo COMPLETION_WAIT_FINISHED
