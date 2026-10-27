#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n probe_environment.sh
# RUN_COMMAND: bash probe_environment.sh
# EXPECTED_SUCCESS_MARKER: ENVIRONMENT_PROBE_COMPLETE
set -euo pipefail
BASE=/mnt/BioDisk/wangmingyuan
printf 'HOST\n'; hostname
printf 'CPU\n'; nproc
printf 'R libraries\n'
find "$BASE/R" "$BASE/project26210_envs" "$BASE/miniconda3/envs/R_4.5.2" -maxdepth 7 -type f \( -path '*/copykat/DESCRIPTION' -o -path '*/infercnv/DESCRIPTION' \) -print 2>/dev/null || true
printf 'Known environment versions\n'
"$BASE/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript" --version
"$BASE/miniconda3/envs/batch_bench/bin/python" -c 'import importlib.util; print({m:bool(importlib.util.find_spec(m)) for m in ["h5py","anndata","scanpy","scipy","matplotlib","pandas"]})'
printf 'Gene coordinate candidates\n'
find "$BASE/ref" "$BASE/project26210_part1_cnv_stage" "$BASE/project26210_envs" -maxdepth 6 -type f \( -iname '*gene*order*' -o -iname '*grch38*' -o -iname '*hg38*' -o -iname '*copykat*' \) -print 2>/dev/null | head -60 || true
echo ENVIRONMENT_PROBE_COMPLETE
