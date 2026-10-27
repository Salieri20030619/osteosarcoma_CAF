#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n run_cnv_cohort.sh
# RUN_COMMAND: bash run_cnv_cohort.sh <infercnv|copykat> <0|1|2>
# EXPECTED_SUCCESS_MARKER: CNV_COHORT_COMPLETE
set -uo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
METHOD=${1:?method}; SHARD=${2:?shard}
case "$METHOD:$SHARD" in
  infercnv:0) samples=(S0059 S0217 S0218 SC069 SC073) ;;
  copykat:1) samples=(S0059 S0218 SC073) ;;
  copykat:2) samples=(S0217 SC069) ;;
  *) exit 2 ;;
esac
echo $$ > "$ROOT/logs/cohort_${METHOD}_${SHARD}.pid"
failures=0
for sample in "${samples[@]}"; do
  if bash "$ROOT/code/run_cnv_job.sh" "$METHOD" "$sample" > "$ROOT/logs/${METHOD}_${sample}.log" 2>&1; then
    echo "SAMPLE_COMPLETE $METHOD $sample"
  else
    echo "SAMPLE_FAILED $METHOD $sample"; failures=$((failures+1))
  fi
done
printf '%s\n' "$failures" > "$ROOT/logs/cohort_${METHOD}_${SHARD}.failures"
if (( failures == 0 )); then echo CNV_COHORT_COMPLETE; else echo CNV_COHORT_FINISHED_WITH_FAILURES; exit 1; fi
