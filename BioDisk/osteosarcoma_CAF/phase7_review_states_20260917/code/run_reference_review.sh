#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n code/run_reference_review.sh
# RUN_COMMAND: bash code/run_reference_review.sh
# EXPECTED_SUCCESS_MARKER: REFERENCE_REVIEW_JOBS_EXITED
set -uo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
exec 8>logs/reference_review.lock
flock -n 8 || exit 73
echo $$ > logs/reference_review.pid
for sample in GSE299023__OS.C.33.BX GSE152048__BC3 GSE152048__BC6; do
 bash code/run_scevan_job.sh "$sample" > "logs/scevan_${sample}.log" 2>&1
 printf '%s\t%s\n' "$sample" "$?" >> logs/reference_review_exit_codes.tsv
done
printf '%s\n' REFERENCE_REVIEW_JOBS_EXITED > REVIEW_JOBS_EXITED
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u code/integrate_identity.py > logs/integrate.log 2>&1
printf '%s\n' "$?" > logs/integrate.exit_code
