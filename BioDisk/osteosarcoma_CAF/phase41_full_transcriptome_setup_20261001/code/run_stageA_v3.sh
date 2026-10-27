#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n run_stageA.sh
# RUN_COMMAND: bash run_stageA.sh
# EXPECTED_SUCCESS_MARKER: STAGE_A_FULL_TRANSCRIPTOME_SETUP_COMPLETE + COMPLETE=SUCCESS + setup.exit_code=0
set -u
out=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase41_full_transcriptome_setup_20261001
echo $$ > "$out/setup_runner.pid"
/usr/bin/time -v /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u "$out/code/stageA_setup_v3.py" > "$out/logs/setup_v3.log" 2>&1
rc=$?
printf '%s\n' "$rc" > "$out/setup_runner.exit_code"
tail -35 "$out/logs/setup_v3.log"
exit "$rc"