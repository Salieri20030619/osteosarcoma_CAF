#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n finish_artifacts.sh
# RUN_COMMAND: bash finish_artifacts.sh
# EXPECTED_SUCCESS_MARKER: PHASE3_ARTIFACTS_COMPLETE
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
SOURCE=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/OS_CAF_phase2
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
RBIN=/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript
echo $$ > "$ROOT/logs/finish_artifacts.pid"
trap 'rc=$?; printf "%s\n" "$rc" > "$ROOT/logs/finish_artifacts.exit_code"' EXIT
test -f "$ROOT/identity/COMPLETE"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 MPLBACKEND=Agg
"$PY" -m py_compile "$ROOT/code/validate_identity.py" "$ROOT/code/build_delivery.py" "$ROOT/code/plot_cnv_summary.py" "$ROOT/code/plot_caf_sets.py"
"$RBIN" --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' "$ROOT/code/export_cnv_plot_data.R"
"$PY" "$ROOT/code/validate_identity.py" --root "$SOURCE" --phase3 "$ROOT" > "$ROOT/logs/identity_validation_final.log" 2>&1
"$RBIN" --vanilla "$ROOT/code/export_cnv_plot_data.R" > "$ROOT/logs/cnv_plot_export.log" 2>&1
"$PY" "$ROOT/code/plot_cnv_summary.py" "$ROOT" > "$ROOT/logs/cnv_plot_summary.log" 2>&1
"$PY" "$ROOT/code/plot_caf_sets.py" "$SOURCE" "$ROOT" > "$ROOT/logs/caf_set_plots.log" 2>&1
"$PY" -m py_compile "$ROOT/code/redraw_marker_plots.py"
"$PY" "$ROOT/code/redraw_marker_plots.py" "$ROOT" > "$ROOT/logs/marker_plot_layout.log" 2>&1
"$PY" -m py_compile "$ROOT/code/write_figure_captions.py"
"$PY" "$ROOT/code/write_figure_captions.py" "$ROOT" > "$ROOT/logs/figure_captions.log" 2>&1
"$PY" "$ROOT/code/build_delivery.py" --root "$SOURCE" --phase3 "$ROOT"
echo PHASE3_ARTIFACTS_COMPLETE
