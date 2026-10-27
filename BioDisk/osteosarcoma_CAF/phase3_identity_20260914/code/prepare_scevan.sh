#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n prepare_scevan.sh
# RUN_COMMAND: bash prepare_scevan.sh
# EXPECTED_SUCCESS_MARKER: SCEVAN_PREP_COMPLETE
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
echo $$ > "$ROOT/logs/scevan_prepare.pid"
trap 'rc=$?; printf "%s\n" "$rc" > "$ROOT/logs/scevan_prepare.exit_code"' EXIT
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
export PATH="/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin:/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin:$PATH"
export MAKEFLAGS=-j4
unset R_MAKEVARS_USER
for item in AntonioDeFalco/SCEVAN:main miccec/yaGST:master; do
  repo=${item%:*}; branch=${item#*:}; name=${repo#*/}
  if [[ ! -s "$ROOT/resources/${name}.tar.gz" ]]; then
    curl -LfsS --retry 2 --max-time 180 "https://api.github.com/repos/$repo/commits/$branch" -o "$ROOT/resources/${name}_commit.json"
    sha=$(/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -c 'import json,sys;print(json.load(open(sys.argv[1]))["sha"])' "$ROOT/resources/${name}_commit.json")
    curl -LfsS --retry 2 --max-time 600 "https://codeload.github.com/$repo/tar.gz/$sha" -o "$ROOT/resources/${name}.tar.gz.part"
    mv "$ROOT/resources/${name}.tar.gz.part" "$ROOT/resources/${name}.tar.gz"
  fi
  sha256sum "$ROOT/resources/${name}.tar.gz"
done
Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' "$ROOT/code/install_scevan.R"
/usr/bin/time -v -o "$ROOT/logs/scevan_prepare_resources.txt" Rscript --vanilla "$ROOT/code/install_scevan.R"
echo SCEVAN_PREP_COMPLETE
