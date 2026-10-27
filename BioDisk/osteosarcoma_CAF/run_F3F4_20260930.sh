#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n run_F3F4_20260930.sh; batch_bench/bin/python -m py_compile review_F3F4_20260930.py
# RUN_COMMAND: bash /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/run_F3F4_20260930.sh
# EXPECTED_SUCCESS_MARKER: F3F4_EXISTING_GENOMEWIDE_REVIEW_COMPLETE; phase40/COMPLETE=SUCCESS; review.exit_code=0
set -u
root=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF
out="$root/phase40_F3F4_matrix_review_20260930"
mkdir -p "$out/logs" "$out/code"
echo $$ > "$out/review.pid"
cp "$root/review_F3F4_20260930.py" "$out/code/review_F3F4.py"
/usr/bin/time -v /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u "$root/review_F3F4_20260930.py" > "$out/logs/review.log" 2>&1
rc=$?
printf '%s\n' "$rc" > "$out/review.exit_code"
tail -95 "$out/logs/review.log"
exit "$rc"