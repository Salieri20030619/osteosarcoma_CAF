#!/usr/bin/env bash
set -uo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
exec 8>remediation_20260917/logs/cnv_reference.lock
flock -n 8 || exit 73
echo $$ > remediation_20260917/logs/cnv_reference.pid
for short in BC10 BC17 BC2 BC5 BC21 BC22; do
 bash remediation_20260917/cnv_immune_reference/code/run_cnv_job.sh infercnv "GSE152048__$short" > "remediation_20260917/logs/cnv_$short.log" 2>&1
 printf '%s\t%s\n' "$short" "$?" >> remediation_20260917/logs/cnv_exit_codes.tsv
done
printf 'CNV_SENSITIVITY_RUNS_EXITED\n' > remediation_20260917/logs/CNV_RUNS_EXITED
