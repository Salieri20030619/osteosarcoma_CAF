#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n switch_to_geo.sh
# RUN_COMMAND: bash switch_to_geo.sh
# EXPECTED_SUCCESS_MARKER: GEO_TRANSPORT_STARTED
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
cd "$ROOT"
bash -n run_geo.sh
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -m py_compile download_validate_geo.py
OLD=$(cat python.pid)
if [[ -r /proc/$OLD/cmdline ]]; then
 CMD=$(tr '\0' ' ' < /proc/$OLD/cmdline)
 if [[ "$CMD" == *"download_validate.py --root $ROOT"* ]]; then
  pkill -TERM -P "$OLD" || true
  kill -TERM "$OLD"
 else
  echo PROCESS_MISMATCH_STOPPED; exit 2
 fi
fi
nohup bash run_geo.sh "$ROOT" > download_geo.log 2>&1 < /dev/null &
echo $! > geo_launch.pid
echo GEO_TRANSPORT_STARTED
