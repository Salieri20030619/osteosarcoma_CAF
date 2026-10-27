#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n start_http1.sh
# RUN_COMMAND: bash start_http1.sh
# EXPECTED_SUCCESS_MARKER: HTTP1_DOWNLOAD_STARTED
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
cd "$ROOT"
bash -n run_http1.sh
bash -n wait_http1.sh
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -m py_compile download_validate_http1.py finalize_download_review.py
nohup bash run_http1.sh "$ROOT" > download_http1.log 2>&1 < /dev/null &
echo $! > http1_launch.pid
echo HTTP1_DOWNLOAD_STARTED
