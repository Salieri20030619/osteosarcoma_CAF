#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n run_http1.sh
# RUN_COMMAND: bash run_http1.sh /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
# EXPECTED_SUCCESS_MARKER: NON_PDX_DOWNLOAD_VALIDATION_COMPLETE
set -euo pipefail
ROOT="$1"
cd "$ROOT"
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$PY" -m py_compile download_validate_http1.py
echo $$ > runner.pid
date -Is > http1_started_at.txt
/usr/bin/time -v -o resource_usage_http1.txt "$PY" -u download_validate_http1.py --root "$ROOT"
date -Is > finished_at.txt
