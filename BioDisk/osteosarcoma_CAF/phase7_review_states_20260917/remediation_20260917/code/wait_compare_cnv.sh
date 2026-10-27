#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
while kill -0 2415198 2>/dev/null; do sleep 900; done
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u remediation_20260917/code/compare_cnv.py > remediation_20260917/logs/cnv_comparison.log 2>&1
