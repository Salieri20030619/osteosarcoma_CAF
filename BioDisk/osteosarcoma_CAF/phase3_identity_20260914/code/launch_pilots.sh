#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n launch_pilots.sh
# RUN_COMMAND: bash launch_pilots.sh
# EXPECTED_SUCCESS_MARKER: CNV_PILOTS_LAUNCHED
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
test -s "$ROOT/resources/coordinate_audit.tsv"
bash -n "$ROOT/code/run_cnv_job.sh"
for sample in SC072 S0058; do
  for method in infercnv copykat; do
    out="$ROOT/$method/$sample"
    if [[ -f "$out/COMPLETE" || -d "$out/job.lock" ]]; then echo "SKIP_EXISTING $method $sample"; continue; fi
    nohup bash "$ROOT/code/run_cnv_job.sh" "$method" "$sample" > "$ROOT/logs/${method}_${sample}.log" 2>&1 < /dev/null &
    printf '%s\t%s\t%s\n' "$method" "$sample" "$!" >> "$ROOT/logs/pilot_launch.tsv"
  done
done
echo CNV_PILOTS_LAUNCHED
