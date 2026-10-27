#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n run_stageB.sh
# RUN_COMMAND: bash run_stageB.sh
# EXPECTED_SUCCESS_MARKER: STAGE_B_F3F4_IDENTITY_STATES_COMPLETE + COMPLETE=SUCCESS + identity_states.exit_code=0
set -u
out=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase42_F3F4_identity_states_20261001
mkdir -p "$out/logs"
echo $$ > "$out/runner.pid"
/usr/bin/time -v /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u "$out/code/stageB_identity_states_v3.py" > "$out/logs/identity_states_v3.log" 2>&1
rc=$?
printf '%s\n' "$rc" > "$out/runner.exit_code"
tail -40 "$out/logs/identity_states_v3.log"
exit "$rc"