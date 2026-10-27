#!/usr/bin/env bash
set -uo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
exec 8>remediation_20260917/logs/decontx.lock
flock -n 8 || exit 73
echo $$ > remediation_20260917/logs/decontx.pid
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate R_4.5.2
for sample in GSE152048__BC10 GSE152048__BC5 GSE152048__BC21 GSE152048__BC2 GSE152048__BC22 GSE152048__BC11 GSE152048__BC17 GSE152048__BC16 GSE152048__BC20 GSE152048__BC3 GSE152048__BC6 GSE299023__001_BX_COM GSE299023__003_BX_COM GSE299023__010_BX_COM GSE299023__OS.C.21.BX.OCT GSE299023__OS.C.32.BX GSE299023__OS.C.33.BX; do
 if [[ ! -f remediation_20260917/decontx_inputs/$sample/COMPLETE ]]; then echo "MISSING_INPUT $sample"; exit 74; fi
 /usr/bin/time -v -o "remediation_20260917/logs/decontx_$sample.resources" flock "remediation_20260917/logs/decontx_$sample.lock" Rscript --vanilla remediation_20260917/code/run_decontx.R "$sample" > "remediation_20260917/logs/decontx_$sample.log" 2>&1
 printf '%s\t%s\n' "$sample" "$?" >> remediation_20260917/logs/decontx_exit_codes.tsv
done
printf 'DECONTX_RUNS_EXITED\n' > remediation_20260917/logs/DECONTX_RUNS_EXITED
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u remediation_20260917/code/compare_decontx.py > remediation_20260917/logs/decontx_comparison.log 2>&1
printf '%s\n' "$?" > remediation_20260917/logs/decontx_comparison.exit_code
