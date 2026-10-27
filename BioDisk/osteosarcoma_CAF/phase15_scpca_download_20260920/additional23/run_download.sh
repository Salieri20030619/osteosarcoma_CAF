#!/usr/bin/env bash
set -u
set -o pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase15_scpca_download_20260920/additional23 || exit 99
py=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
export SCPCA_DOWNLOAD_PROXY=http://127.0.0.1:7897
"$py" -m py_compile range_zip.py download_parallel.py || exit 98
"$py" -u download_parallel.py 2>&1 | tee download_proxy.log
code=$?
printf '%s\n' "$code" > exit_code
if [ "$code" -ne 0 ]; then exit "$code"; fi
"$py" verify_download.py > verification.log 2>&1
qa=$?
printf '%s\n' "$qa" > qa_exit_code
exit "$qa"
