#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n wait_completion.sh
# RUN_COMMAND: bash wait_completion.sh
# EXPECTED_SUCCESS_MARKER: NON_PDX_DOWNLOAD_VALIDATION_COMPLETE
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
PID=$(cat "$ROOT/geo_launch.pid")
while kill -0 "$PID" 2>/dev/null; do
 if [[ -r /proc/$PID/stat ]] && [[ $(awk '{print $3}' /proc/$PID/stat) == Z ]]; then break; fi
 sleep 900
done
tail -n 25 "$ROOT/download_geo.log"
if [[ -f "$ROOT/COMPLETE" ]]; then cat "$ROOT/COMPLETE"; else echo DOWNLOAD_OR_VALIDATION_FAILED; exit 1; fi
