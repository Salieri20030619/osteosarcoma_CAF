#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n prepare_resources_v3.sh
# RUN_COMMAND: bash prepare_resources_v3.sh
# EXPECTED_SUCCESS_MARKER: RESOURCES_READY
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
source /mnt/BioDisk/wangmingyuan/miniconda3/etc/profile.d/conda.sh
conda activate /mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2
unset R_MAKEVARS_USER
export PATH="/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin:/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin:$PATH"
export MAKEFLAGS=-j4
echo $$ > "$ROOT/logs/resources_v3.pid"
trap 'rc=$?; if (( rc != 0 )); then printf "%s\n" "$rc" > "$ROOT/logs/RESOURCES_V3_FAILED"; fi' EXIT
command -v x86_64-conda-linux-gnu-gcc x86_64-conda-linux-gnu-c++ x86_64-conda-linux-gnu-gfortran
x86_64-conda-linux-gnu-gcc --version | head -1
fetch() {
 local url="$1" name="$2"
 if [[ ! -s "$ROOT/resources/$name" ]]; then
   curl -L --fail --retry 2 --connect-timeout 20 --max-time 900 "$url" -o "$ROOT/resources/$name.part"
   mv "$ROOT/resources/$name.part" "$ROOT/resources/$name"
 fi
 sha256sum "$ROOT/resources/$name"
}
fetch https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz Homo_sapiens.GRCh38.84.gtf.gz &
P38=$!
fetch https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz Homo_sapiens.GRCh37.75.gtf.gz &
P37=$!
Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' "$ROOT/code/install_copykat.R"
Rscript --vanilla "$ROOT/code/install_copykat.R"
wait "$P38"; wait "$P37"
gzip -t "$ROOT/resources/Homo_sapiens.GRCh38.84.gtf.gz" "$ROOT/resources/Homo_sapiens.GRCh37.75.gtf.gz"
date --iso-8601=seconds > "$ROOT/logs/RESOURCES_READY"
echo RESOURCES_READY
